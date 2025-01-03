version 1.1

task active_regions_merge {
    meta {
        description: "Merge active regions"
        outputs: {
            combined_promoters: "Merged active regions"
        }
    }

    parameter_meta {
        promoters: "Promoter regions in BED format"
        loop_bed: "Loop anchors in BED format"
        bam: "BAM file"
    }

    input {
        File promoters
        File loop_bed
        File bam
    }

    String outfile = basename(promoters, ".bed") + ".LoopAnchors.Enhs.combined.sort.bed"

    #@ except: LineWidth
    command <<<
        prom=~{basename(promoters,".bed")}
        sort \
            -k1,1 \
            -k2,2n \
            -k4,4 \
            ~{promoters} | \
            awk \
            -F'\t' \
            -v OFS="\t" \
            '{ a[$1"\t"$2"\t"$3]=($1"\t"$2"\t"$3 in a? a[$1"\t"$2"\t"$3]",":"")$4 }END{ for(i in a) print i,a[i] }' | \
            sort \
            -k1,1 \
            -k2,2n \
            > $prom".unique.bed"

        #get read counts
        intersectBed -c -a $prom".unique.bed" -b ~{bam} > $prom"_Signal"

        #compute top 2/3
        length=$(wc -l $prom"_Signal" | awk -F\  '{print $1}')
        keepN=$((length * 2 / 3))

        sort -k5,5nr $prom"_Signal" \
            | head -n "$keepN" \
            | sort -k1,1 -k2,2n \
            | cut -f 1-3 \
            > "${prom}_Signal.keep.${keepN}.3cols.bed"

        #get non-promoter region
        cut -f 1-3 ~{loop_bed} | \
            sort -k1,1 -k2,2n \
            > ~{basename(loop_bed)}
        bedtools \
            subtract \
            -a ~{basename(loop_bed)} \
            -b $prom".unique.bed" | \
            sort \
            > Enhancers.Subtract.4kbPromoters.sort.bed

        cat \
            "${prom}_Signal.keep.${keepN}.3cols.bed" \
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
        container: "ghcr.io/stjudecloud/bedtools:2.31.1"
    }
}

task extract_promoters {
    meta {
        description: "Extract promoters from a GTF file"
        outputs: {
            promoter: "Promoter regions in BED format"
        }
    }

       parameter_meta {
        annotation: "GTF (optionally gzip compressed) file containing gene annotations"
    }

    input {
        File annotation
    }

    String base = basename(annotation, ".gz")

    #@ except: LineWidth
    command <<<
        set -euo pipefail

        gunzip -c ~{annotation} > ~{base} \
           || ln -sf ~{annotation} ~{base}


        sed -i 's/ /\t/g' ~{base}

        # Find all transcripts, excluding the mitochondrial chromosome
        # Add a 2000bp buffer to the transcription start site to include promoter
        # Extract gene_id and transcript_id ($10 and $12).
        # Subtract the suffix from the trasncript ID.
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

    #@ except: LineWidth
    command <<<
        set -euo pipefail

        gunzip -c ~{annotation} > ~{base} \
           || ln -sf ~{annotation} ~{base}

        sed -i 's/ /\t/g' ~{base}

        # Find all transcripts, excluding the mitochondrial chromosome
        # Extract gene_id and transcript_id ($10 and $12).
        # Subtract the suffix from the trasncript ID.
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
        mapping_stats_read1: "HiC-Pro mapping statistics for read 1"
        mapping_stats_read2: "HiC-Pro mapping statistics for read 2"
        pairing_stats: "HiC-Pro pairing statistics"
        peaks_bed: "Peaks in BED format"
        fithichip_bed: "FithiChIP interactions in BED format"
        fithichip_q01_bed: "FithiChIP interactions with q-value < 0.01 in BED format"
        prefix: "Prefix for the output file"
    }

    input {
        File all_valid_pairs_stats
        File mapping_stats_read1
        File mapping_stats_read2
        File pairing_stats
        File? peaks_bed
        File? fithichip_bed
        File? fithichip_q01_bed
        String prefix = basename(all_valid_pairs_stats, "_allValidPairs.mergestat")
    }

    Int disk_size_gb = ceil(
        size(all_valid_pairs_stats, "GiB") +
        size(mapping_stats_read1, "GiB") +
        size(mapping_stats_read2, "GiB") +
        size(pairing_stats, "GiB") +
        size(peaks_bed, "GiB") +
        size(fithichip_bed, "GiB") +
        size(fithichip_q01_bed, "GiB")
    ) + 2

    command <<<
        python /usr/local/bin/qc_hic.py \
            --all_valid_pairs_stats ~{all_valid_pairs_stats} \
            --mapping_stats_read1 ~{mapping_stats_read1} \
            --mapping_stats_read2 ~{mapping_stats_read2} \
            --pairing_stats ~{pairing_stats} \
            ~{if defined(peaks_bed) then "--peaks_bed ~{peaks_bed}" else ""} \
            ~{if defined(fithichip_bed) then "--fithichip_bed ~{fithichip_bed}" else ""} \
            ~{(
                if defined(fithichip_q01_bed)
                then "--fithichip_q01_bed ~{fithichip_q01_bed}"
                else ""
            )} \
            --prefix ~{prefix}
    >>>

    output {
        File qc_report = prefix + "_QCreport.txt"
    }

    runtime {
        container: "ghcr.io/stjudecloud/hilow:1.0.0"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        memory: "4 GB"
        maxRetries: 1
    }
}

task filter {
    meta {
        description: "Filter Hi-C results using an exclude list"
        outputs: {
            filtered_pairs: "Filtered pairs file",
            removed_pairs: "Removed pairs file",
            filtered_stats: "Filtered statistics file",
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
    Int disk_size_gb = ceil(size(all_valid_pairs, "GiB")) + 2

    #@ except: LineWidth
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

        all=$(wc -l ~{all_valid_pairs} |cut -d " " -f 1)
        filtered=$(wc -l ~{prefix}.allValidPairs.removed | cut -d " " -f 1)
        percentage=$(echo "$filtered/$all*100"|bc -l)
        echo -e "$filtered read pairs are filtered out from a total of $all, accounting for ${percentage} %\n" \
            | tee \
            > ~{prefix}.allValidPairs.stats

    >>>

    output {
        File filtered_pairs = prefix + ".allValidPairs.filtered"
        File removed_pairs = prefix + ".allValidPairs.removed"
        File filtered_stats = prefix + ".allValidPairs.stats"
    }

    runtime {
        container: "ghcr.io/stjudecloud/bedtools:branch-hic_workflow-2.31.1"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        memory: "4 GB"
        maxRetries: 1
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
    Int disk_size_gb = ceil(size(all_valid_pairs, "GiB")) * 2

    command <<<
        set -euo pipefail

        /HiC-Pro_3.0.0/bin/utils/hicpro2juicebox.sh \
            -i ~{all_valid_pairs} \
            -g ~{chromsizes} \
            -j /opt/juicer-1.6.2/CPU/common/juicer_tools.1.7.6_jcuda.0.8.jar
    >>>

    output {
        File hic_file = prefix + ".hic"
    }

    runtime {
        container: "ghcr.io/stjudecloud/juicertools:1.6.2"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        memory: "18 GB"
        maxRetries: 1
    }
}
