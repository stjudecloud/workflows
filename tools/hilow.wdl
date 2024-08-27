version 1.1

task activeRegionsMerge {
    meta {}
    parameter_meta {}

    input {
        File promoters
        File loopBed
        File bam
        String prefix
    }

    String outfile = basename(promoters, ".bed") + ".LoopAnchors.Enhs.combined.sort.bed"

    command <<<
        prom=~{basename(promoters,".bed")}
        sort \
            -k1,1 \
            -k2,2n \
            -k4,4 \
            ~{promoters} \
            | awk \
            -F'\t' \
            -v OFS="\t" \
            '{ a[$1"\t"$2"\t"$3]=($1"\t"$2"\t"$3 in a? a[$1"\t"$2"\t"$3]",":"")$4 }END{ for(i in a) print i,a[i] }'
            | sort \
            -k1,1 \
            -k2,2n \
            > $prom".unique.bed"

        #get read counts
        intersectBed -c -a $prom".unique.bed" -b ~{bam} > $prom"_Signal"

        #compute top 2/3
        length=$(wc -l $prom"_Signal" | awk -F\  '{print $1}'); 
        keepN=$(echo "$(($length * 2 / 3))"); 

        sort -k5,5nr $prom"_Signal" \
            | head -n $keepN \
            | sort -k1,1 -k2,2n \
            | cut -f 1-3 \
            > $prom"_Signal.keep."$keepN".3cols.bed"

        #get non-promoter region
        cut -f 1-3 ~{loopBed} 
            | sort -k1,1 -k2,2n \
            > ~{basename(loopBed)}
        bedtools \
            subtract \
            -a ~{basename(loopBed)} \
            -b $prom".unique.bed" \
            | sort \
            > Enhancers.Subtract.4kbPromoters.sort.bed

        cat \
            $prom"_Signal.keep."$keepN".3cols.bed" \
            Enhancers.Subtract.4kbPromoters.sort.bed \
            | sort -k1,1 -k2,2n \
            > ~{outfile}
    >>>

    output {
        File combined_promoters = outfile
    }

    runtime {
        cpu: 8
        memory: "20 GB"
        maxRetries: 1
        container: "adthrasher/bedtools:2.31.1"
    }
}

task extract_promoters {
    meta {
        description: "Extract promoters from a GTF file"
        outputs: {
            promoters: "Promoter regions in BED format"
        }
    }

       parameter_meta {
        annotation: "GTF (optionally gzip compressed) file containing gene annotations"
    }

    input {
        File annotation
    }

    String base = basename(annotation, ".gz")

    command <<<
        set -euo pipefail

        gunzip -c ~{annotation} > ~{base} \
           || ln -sf ~{annotation} ~{base}


        sed -i 's/ /\t/g' ~{base}

        # Find all transcripts, excluding the mitochondrial chromosome
        # Add a 2000bp buffer to the transcription start site to include promoter
        # Extract gene_id and transcript_id ($10 and $12). Subtract the suffix from the trasncript ID.
        # Convert to BED format
        awk -F\\t '{ if($3 == "transcript" && $1 !~ "chrM")
            if($7 =="+")
                print $1 "++" $4-2000 "++" $4+2000 "++" substr($12,2,length($12)-5) "|" substr($10,2,length($10)-3)
            else if($7=="-")
                print $1 "++" $5-2000 "++" $5+2000 "++" substr($12,2,length($12)-5) "|" substr($10,2,length($10)-3) }' \
            ~{base} \
            | grep -v "unassigned" \
            | sort \
            | uniq \
            | sed s/\+\+/\\t/g \
            | sort -k1,1 -k2,2n > promoter_regions.bed
    >>>

    output {
        File promoter = "promoter_regions.bed"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
        container: "ubuntu:22.04"
    }
}

task extract_genes {
    meta {
        description: "Extract genes from a GTF file"
        outputs: {
            genes: "Gene regions in BED format"
        }
    }

    parameter_meta {
        annotation: "GTF (optionally gzip compressed) file containing gene annotations"
    }

    input {
        File annotation
    }

    String base = basename(annotation, ".gz")

    command <<<
        set -euo pipefail

        gunzip -c ~{annotation} > ~{base} \
           || ln -sf ~{annotation} ~{base}

        sed -i 's/ /\t/g' ~{base}

        # Find all transcripts, excluding the mitochondrial chromosome
        # Extract gene_id and transcript_id ($10 and $12). Subtract the suffix from the trasncript ID.
        # Convert to BED format
        awk -F\\t '{ if($3 == "transcript" && $1 !~ "chrM") 
            print $1 "++" $4 "++" $5 "++" substr($12,2,length($12)-5) "|" substr($10,2,length($10)-3)}' \
            ~{base} \
            | grep -v "unassigned" \
            | sort \
            | uniq \
            | sed s/\+\+/\\t/g \
            | sort -k1,1 -k2,2n > gene_regions.bed
    >>>

    output {
        File genes = "gene_regions.bed"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
        container: "ubuntu:22.04"
    }
}


task qcreport {
    meta {
        description: "Write QC report for Hi-C experiment"
        outputs: {
            qc_report: "QC report file"
        }
    }

    parameter_meta {
        all_valid_pairs_stats: "HiC-Pro allValidPairs mergestat statistics file"
        mapping_stats_R1: "HiC-Pro mapping statistics for read 1"
        mapping_stats_R2: "HiC-Pro mapping statistics for read 2"
        pairing_stats: "HiC-Pro pairing statistics"
        peaks_bed: "Peaks in BED format"
        fithichip_bed: "FithiChIP interactions in BED format"
        fithichip_q01_bed: "FithiChIP interactions with q-value < 0.01 in BED format"
        prefix: "Prefix for the output file"
    }

    input {
        File all_valid_pairs_stats
        File mapping_stats_R1
        File mapping_stats_R2
        File pairing_stats
        File? peaks_bed
        File? fithichip_bed
        File? fithichip_q01_bed
        String prefix = basename(all_valid_pairs_stats, "_allValidPairs.mergestat")
    }

    command <<<
        python <<CODE
        import os
        FILES = ["~{all_valid_pairs_stats}", "~{pairing_stats}", "~{mapping_stats_R1}", "~{mapping_stats_R2}"]
        VARIABLE = ["valid_interaction", "cis_shortRange", "cis_longRange", "Total_pairs_processed", "Reported_pairs", "total_R2", "mapped_R2", "total_R1", "mapped_R1"]
        RESULTS = {}
        for eachfile in FILES:
            with open(eachfile) as f:
                for line in f:
                    array = line.split()
                    if array[0] in VARIABLE:
                        RESULTS[array[0]] = int(array[1])

        PERCENTAGES = {}
        COUNTS = {}
        VERDICT = {}
        REP_VAR = ["R1_aligned", "R2_aligned", "valid_interactionPairs", "cis_shortRange", "cis_longRange"]

        if os.path.isfile('~{fithichip_bed}'):
            with open('~{fithichip_q01_bed}') as fithichip:
                LOOPS_SIGNIFICANT = len(fithichip.readlines())-1
            with open('~{fithichip_bed}') as fithichip:
                LOOPS = len(fithichip.readlines())-1
        if os.path.isfile("~{peaks_bed}"):
            with open("~{peaks_bed}") as peaks:
                PEAKS = len(peaks.readlines())-1

        PERCENTAGES["R1_aligned"] = round((RESULTS["mapped_R1"]*100)/RESULTS["total_R1"])
        PERCENTAGES["R2_aligned"] = round((RESULTS["mapped_R2"]*100)/RESULTS["total_R2"])
        PERCENTAGES["valid_interactionPairs"] = round((RESULTS["valid_interaction"]*100)/RESULTS["Reported_pairs"])
        PERCENTAGES["cis_shortRange"] = round((RESULTS["cis_shortRange"]*100)/RESULTS["valid_interaction"])
        PERCENTAGES["cis_longRange"] = round((RESULTS["cis_longRange"]*100)/RESULTS["valid_interaction"])
        COUNTS["R1_aligned"] = RESULTS["mapped_R1"]
        COUNTS["R2_aligned"] = RESULTS["mapped_R2"]
        COUNTS["valid_interactionPairs"] = RESULTS["valid_interaction"]
        COUNTS["cis_shortRange"] = RESULTS["cis_shortRange"]
        COUNTS["cis_longRange"] = RESULTS["cis_longRange"]

        for aligned in ["R1_aligned", "R2_aligned"]:
            if PERCENTAGES[aligned] > 80:
                VERDICT[aligned] = "GOOD"
            else:
                VERDICT[aligned] = "BAD"

        if PERCENTAGES["valid_interactionPairs"] > 50:
            VERDICT["valid_interactionPairs"] = "GOOD"
        else:
            VERDICT["valid_interactionPairs"] = "BAD"

        if PERCENTAGES["cis_shortRange"] > 50:
            VERDICT["cis_shortRange"] = "BAD"
        elif PERCENTAGES["cis_shortRange"] > 30:
            VERDICT["cis_shortRange"] = "MARGINAL"
        else:
            VERDICT["cis_shortRange"] = "GOOD"

        if PERCENTAGES["cis_longRange"] > 40:
            VERDICT["cis_longRange"] = "GOOD"
        elif PERCENTAGES["cis_longRange"] > 20:
            VERDICT["cis_longRange"] = "MARGINAL"
        else:
            VERDICT["cis_longRange"] = "BAD"

        REPORT = open("~{prefix}_QCreport.txt", 'w')
        REPORT.write("STAT\tCOUNTS\tPERCENTAGE\tVERDICT\n")
        REPORT.write("Total_pairs_processed\t" + str(RESULTS["Total_pairs_processed"]) + "\n")
        for var in REP_VAR:
            REPORT.write(var + "\t" + str(COUNTS[var]) + "\t" + str(PERCENTAGES[var]) + "\t" + VERDICT[var] + '\n')

        if os.path.isfile("~{peaks_bed}"):
            REPORT.write("peaks\t" + str(PEAKS) + "\n")
        if os.path.isfile('~{fithichip_bed}'):
            REPORT.write("loops\t" + str(LOOPS) + "\n")
            REPORT.write("loops_significant\t" + str(LOOPS_SIGNIFICANT) + "\n")
        CODE
    >>>

    output {
        File qc_report = prefix + "_QCreport.txt"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
        container: "python:3.9.19"
    }
}


task filter {
    meta {
        description: "Filter Hi-C results using an exclude list"
        outputs: {
            filtered_pairs: "Filtered pairs file",
            filtered_stats: "Filtered statistics file",
            removed_pairs: "Removed pairs file",
        }
    }

    parameter_meta {
        all_valid_pairs: "Valid pairs file"
        chromsizes: "Chromosome sizes file"
        exclude_list: "List of pairs to exclude"
        padding: "Number of bases to add to each feature"
    }

    input {
        File all_valid_pairs
        File chromsizes
        File exclude_list
        Int padding = 50

    }

    String base = basename(exclude_list, ".gz")
    String prefix = basename(all_valid_pairs, ".AllValidPairs")

    command <<<
        set -euo pipefail

        gunzip -c ~{exclude_list} > ~{base} \
           || ln -sf ~{exclude_list} ~{base}


        #left bed
        awk -F "\t" -v OFS="\t" '{print $2,$3,$3,$1}' ~{all_valid_pairs} \
        | slopBed -i stdin -b ~{padding} -g ~{chromsizes} \
        | intersectBed -a stdin -b ~{base} -u > left.bed

        #right bed
        awk -F "\t" -v OFS="\t" '{print $5,$6,$6,$1}' ~{all_valid_pairs} \
        | slopBed -i stdin -b ~{padding} -g ~{chromsizes} \
        | intersectBed -a stdin -b ~{base} -u > right.bed

        cat <(cut -f 4 left.bed) <(cut -f 4 right.bed)|sort -u > filter.pair

        python <<CODE
            from collections import defaultdict
            f=open("filter.pair")
            blackID=defaultdict(int)
            while True:
                line=f.readline()
                if not line:
                    break
                cols=line.strip().split("\t")
                ID=cols[0]
                blackID[ID]=1

            f2=open("~{all_valid_pairs}")
            outfile1="~{prefix}.allValidPairs.filtered"
            outfile2="~{prefix}.allValidPairs.removed"
            of1=open(outfile1,'w')
            of2=open(outfile2,'w')
            while True:
                line=f2.readline()
                if not line:
                    break
                cols=line.strip().split("\t")
                id=cols[0]
                if id in blackID:
                    of2.write(line)
                else:
                    of1.write(line)
        CODE

        all=`wc -l ~{all_valid_pairs} |cut -d " " -f 1`
        filtered=`wc -l ~{prefix}.allValidPairs.removed | cut -d " " -f 1`
        percentage=`echo "$filtered/$all*100"|bc -l`
        echo -e "$filtered read pairs are filtered out from a total of $all, accounting for ${percentage} %\n" | tee > ~{prefix}.allValidPairs.stats

    >>>

    output {
        File filtered_pairs = prefix + ".allValidPairs.filtered"
        File removed_pairs = prefix + ".allValidPairs.removed"
        File filtered_stats = prefix + ".allValidPairs.stats"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
        container: "adthrasher/bedtools:2.31.1"
    }
}

task converthic {
    meta {
        description: "Convert HiC-Pro output (allValidPairs file) to .hic format file"
        outputs: {
            hic_file: "Hi-C file in .hic format"
        }
    }

    parameter_meta {
        all_valid_pairs: "Valid pairs file"
        chromsizes: "Chromosome sizes file"
    }

    input {
        File all_valid_pairs
        File chromsizes
    }

    String prefix = basename(all_valid_pairs)

    command <<<
        set -euo pipefail

        /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh -i ~{all_valid_pairs} -g ~{chromsizes} -j /opt/juicer-1.6.2/CPU/common/juicer_tools.1.7.6_jcuda.0.8.jar
    >>>

    output {
        File hic_file = prefix + ".hic"
    }

    runtime {
        cpu: 1
        memory: "18 GB"
        maxRetries: 1
        container: "adthrasher/juicertools:1.6.2"
    }
}

task qc_hic {
    meta {
        description: "Plot Hi-C quality control statistics"
        outputs: {
            mapping_stats_plot: "Plot of mapping statistcs",
            pairing_stats_plot: "Plot of fragment pairing",
            filtering_stats_plot: "Plot of Hi-C fragments",
            filtering_size_plot: "Plot of Hi-C fragment sizes",
            contacts_stats_plot: "Plot of Hi-C contact ranges",
        }
    }

    parameter_meta {
        mapping_stats: "Mapping statistics files"
        pairing_stats: "Pairing statistics files"
        fragment_stats: "Fragment statistics files"
        contacts_stats: "Contacts statistics files"
        plot_type: {
            description: "Type of plot to generate",
            choices: [
                "all",
                "mapping",
                "pairing",
                "filtering",
                "contacts"
            ]
        }
    }

    input {
        Array[File] mapping_stats
        Array[File] pairing_stats
        Array[File] fragment_stats
        Array[File] contacts_stats
        String plot_type = "all"
        String sample_name = "sample"
        Boolean remove_singleton = false
        Boolean remove_multimapper = false
    }

    Int rmSingle_arg = if remove_singleton then 1 else 0
    Int rmMulti_arg = if remove_multimapper then 1 else 0

    command <<<
        set -euo pipefail

        mkdir plots

        mkdir bowtie
        ln -sf ~{sep(" ", mapping_stats)} bowtie/
        ln -sf ~{sep(" ", pairing_stats)} bowtie/

        mkdir hic
        ln -sf ~{sep(" ", fragment_stats)} hic/

        mkdir stats
        ln -sf ~{sep(" ", contacts_stats)} stats/

        # mapping plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "mapping" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' bwtDir='bowtie' sampleName='~{sample_name}' r1tag='.R1' r2tag='.R2'" \
                /HiC-Pro_3.0.0/scripts/plot_mapping_portion.R \
                plot_mapping_portion.Rout
        fi

        # pairing plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "pairing" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' bwtDir='bowtie' sampleName='~{sample_name}' rmMulti='~{rmMulti_arg}' rmSingle='~{rmSingle_arg}'" \
                /HiC-Pro_3.0.0/scripts/plot_pairing_portion.R \
                plot_pairing_portion.Rout
        fi
        
        # filtering plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "filtering" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' hicDir='hic' sampleName='~{sample_name}' rmMulti='~{rmMulti_arg}' rmSingle='~{rmSingle_arg}'" \
                /HiC-Pro_3.0.0/scripts/plot_hic_fragment.R \
                plot_hic_fragment.Rout
        fi

        # contacts plots
        if [[ "~{plot_type}" == "all" || "~{plot_type}" == "contacts" ]]
        then
            R CMD BATCH \
                --no-save \
                --no-restore \
                "--args picDir='plots' hicDir='hic' statsDir='stats' sampleName='~{sample_name}'" \
                /HiC-Pro_3.0.0/scripts/plot_hic_contacts.R \
                plot_hic_contacts.Rout
        fi
    >>>

    output {
        File? mapping_stats_plot = glob("plots/plotMapping_*.pdf")[0]
        File? pairing_stats_plot = glob("plots/plotMappingPairing_*.pdf")[0]
        File? filtering_stats_plot = glob("plots/plotHiCFragment_*.pdf")[0]
        File? filtering_size_plot = glob("plots/plotHiCFragmentSize_*.pdf")[0]
        File? contacts_stats_plot = glob("plots/plotHiCContactRanges_*.pdf")[0]
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
    }
}
