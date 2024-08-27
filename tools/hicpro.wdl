version 1.1

task cutsite_trimming {
    meta {
        description: "Trim cutsite from reads"
        outputs: {
            cutsite_trimmed_fastq: "FASTQ file with cutsites trimmed",
            cutsite_trimmed_log: "Log file of the trimming process",
        }
    }

    parameter_meta {
        fastq: "Input FASTQ file"
        cutsite: "Cutsite to trim"
    }

    input {
        File fastq
        String cutsite
    }

    String prefix = sub(basename(fastq), ".fq.gz|.fastq.gz|.fq|.fastq", "")

    command <<<
        set -euo pipefail

        fastq=~{basename(fastq, ".gz")}
        gunzip -c ~{fastq} > "$fastq" \
            || ln -sf ~{fastq} "$fastq"

        /HiC-Pro_3.0.0/scripts/cutsite_trimming \
            --fastq $fastq \
            --cutsite ~{cutsite} \
            --out ~{prefix}_cutsite_trimmed.fastq \
            > ~{prefix}_readsTrimming.log \
            2>&1

        gzip ~{prefix}_cutsite_trimmed.fastq && rm $fastq
    >>>

    output {
        File cutsite_trimmed_fastq = prefix + "_cutsite_trimmed.fastq.gz"
        File cutsite_trimmed_log = prefix + "_readsTrimming.log"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
    }
}

task bowtie_pairing {
    meta {
        description: "Merge R1 and R2 BAM files into one paired-end BAM file"
        outputs: {
            combined_bam: "Combined BAM file",
            allele_specific_bam: "Allele specific BAM file",
            combined_stats: "Paired-end statistics file",
        }
    }

    parameter_meta {
        read1_bam: "BAM file containing read 1"
        read2_bam: "BAM file containing read 2"
        allele_specific_tag: "Tag to add to the output"
        prefix: "Prefix for the output file"
        remove_singleton: "Remove singleton reads"
        remove_multimapper: "Remove multi-mapped reads"
    }

    input {
        File read1_bam
        File read2_bam
        File? allele_specific_tag
        String prefix = basename(read1_bam, ".bwt2merged.bam")
        Boolean remove_singleton = false
        Boolean remove_multimapper = false
        Int min_mapq = 10
    }

    command <<<
        set -euo pipefail
        # merge pairs
        python /HiC-Pro_3.0.0/scripts/mergeSAM.py \
            -q ~{min_mapq} \
            -t \
            -v \
            ~{if remove_singleton then "" else "-s"} \
            ~{if remove_multimapper then "" else "-m"} \
            -f ~{read1_bam} \
            -r ~{read2_bam} \
            -o ~{prefix}.bwt2pairs.bam

        ~{if defined(allele_specific_tag)
            then "python /HiC-Pro_3.0.0/scripts/addTagToBAM.py -s ~{allele_specific_tag} -v -r -i ~{prefix}.bwt2pairs.bam -o ~{prefix}.bwt2pairs_allspe.bam"
            else ""
        }

    >>>

    output {
        File combined_bam = prefix + ".bwt2pairs.bam"
        File? allele_specific_bam = prefix + ".bwt2pairs_allspe.bam"
        File combined_stats = prefix + ".bwt2pairs.pairstat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task mapping_stats {
    meta {
        description: "Compute mapping statistics"
        outputs: {
            mapping_stats: "Summarized mapping statistics file"
        }
    }

    parameter_meta {
        combined_bam: "Combined BAM file"
        global_bam: "End-to-end alignments BAM file"
        local_bam: "Local alignments BAM file"
        tag: "Tag to add to the output"
        prefix: "Prefix for the output file"
    }

    input {
        File combined_bam
        File global_bam
        File? local_bam
        String? tag
        String prefix = basename(combined_bam, ".bwt2glob.bam")
    }

    command <<<
        set -euo pipefail

        echo "## HiC-Pro Mapping Statistics" > ~{prefix}.mapstat
        echo "## ~{prefix}.mapstat" >> ~{prefix}.mapstat

        total_reads=$(samtools view -c ~{combined_bam})

        mapped_reads=$(samtools view -c -F 4 ~{combined_bam})

        global_mapped_reads=$(samtools view -c -F 4 ~{global_bam})

        if ~{if defined(local_bam) then true else false}
        then
            local_mapped_reads=$(samtools view -c -F 4 ~{local_bam})
        fi

        echo -e "total_~{tag}\t$total_reads" >> ~{prefix}.mapstat
        echo -e "mapped_~{tag}\t$mapped_reads" >> ~{prefix}.mapstat
        echo -e "global_~{tag}\t$global_mapped_reads" >> ~{prefix}.mapstat
        echo -e "local_~{tag}\t$local_mapped_reads" >> ~{prefix}.mapstat
    >>>

    output {
        File mapping_stats = prefix + ".mapstat"
    }

    runtime {
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task mapped_2hic_fragments {
    meta {
        description: "Keep only valid 3C products"
        outputs: {
            valid_pairs: "Valid pairs file",
            rs_stats: "RS statistics file",
        }
    }

    parameter_meta {
        fragment: "Restriction fragment file in GFF3 format"
        mapped_reads: "Mapped reads in BAM/SAM format"
        shortest_insert_size: {
            description: "Shortest insert size of mapped reads to consider",
            hicpro_field: "MIN_INSERT_SIZE",
        }
        longest_insert_size: {
            description: "Longest insert size of mapped reads to consider",
            hicpro_field: "MAX_INSERT_SIZE",
        }
        shortest_fragment_length: {
            description: "Shortest restriction fragment length to consider",
            hicpro_field: "MIN_FRAG_SIZE",
        }
        longest_fragment_length: {
            description: "Longest restriction fragment length to consider",
            hicpro_field: "MAX_FRAG_SIZE",
        }
        min_cis_distance: {
            description: "Minimum distance between intrachromosomal contact to consider. Filters contacts below this distance. Mainly useful for DNase Hi-C",
            hicpro_field: "MIN_CIS_DIST",
        }
        genotype_tag: "Adds XA tag to report in the valid pairs output for allele specific classification"
        addl_output: {
            description: "Write all additional output files, with information about the discarded reads (self-circle, dangling end, etc.)",
            hicpro_field: "GET_ALL_INTERACTION_CLASSES",
        }
        sam: {
            description: "Output an additional SAM file with flag 'CT' for pairs classification",
            hicpro_field: "GET_PROCESS_BAM",
        }
    }

    input {
        File mapped_reads
        File? fragment
        Boolean genotype_tag = false
        Boolean addl_output = false
        Boolean sam = false
        Int shortest_insert_size = 0
        Int longest_insert_size = 0
        Int shortest_fragment_length = 0
        Int longest_fragment_length = 0
        Int min_cis_distance = 0
    }

    String outfile_name = basename(mapped_reads, ".bam") + ".validPairs"

    command <<<
        set -euo pipefail

        if ~{if defined(fragment) then true else false}
        then
            frag=~{basename(select_first([fragment, ""]), ".gz")}
            gunzip -c ~{fragment} > "$frag" \
                || ln -sf ~{fragment} "$frag"
        fi
        # opts logic from HiC-Pro
        # opts="-v"
        # if [[ $mode == "RS" ]]; then
        #     if [[ "${GET_PROCESS_SAM}" -eq "1" ]]; then opts=$opts" -S"; fi
        #     if [[ "${MIN_FRAG_SIZE}" -ge "0" && "${MIN_FRAG_SIZE}" -ne "" ]]; then opts=$opts" -t ${MIN_FRAG_SIZE}"; fi
        #     if [[ "${MAX_FRAG_SIZE}" -ge "0" && "${MAX_FRAG_SIZE}" -ne "" ]]; then opts=$opts" -m ${MAX_FRAG_SIZE}"; fi
        #     if [[ "${MIN_INSERT_SIZE}" -ge "0" && "${MIN_INSERT_SIZE}" -ne "" ]]; then opts=$opts" -s ${MIN_INSERT_SIZE}"; fi
        #     if [[ "${MAX_INSERT_SIZE}" -ge "0" && "${MAX_INSERT_SIZE}" -ne "" ]]; then opts=$opts" -l ${MAX_INSERT_SIZE}"; fi
        # fi
        # if [[ "${GET_ALL_INTERACTION_CLASSES}" -eq "1" ]]; then opts=$opts" -a"; fi
        # if [[ "${MIN_CIS_DIST}" -ge "0" && "${MIN_CIS_DIST}" -ne "" ]]; then opts=$opts" -d ${MIN_CIS_DIST}"; fi
        # if [[ ! -z ${ALLELE_SPECIFIC_SNP} ]]; then opts=$opts" -g XA"; fi
        ~{if defined(fragment)
            then "python /HiC-Pro_3.0.0/scripts/mapped_2hic_fragments.py -f $frag \\"  # Fragment file found
            else "python /HiC-Pro_3.0.0/scripts/mapped_2hic_dnase.py \\"  # DNAse
        }
        -v \
        -r ~{mapped_reads} \
        ~{if defined(fragment) && sam then "-S" else ""} \
        ~{if addl_output then "-a" else ""} \
        ~{if min_cis_distance > 0 then "-d " + min_cis_distance else ""} \
        ~{if genotype_tag then "-g XA" else ""} \
        ~{if defined(fragment) && shortest_insert_size > 0
            then "-s " + shortest_insert_size
            else ""
        } \
        ~{if defined(fragment) && longest_insert_size > 0
            then "-l " + longest_insert_size
            else ""
        } \
        ~{if defined(fragment) && shortest_fragment_length > 0
            then "-t " + shortest_fragment_length
            else ""
        } \
        ~{if defined(fragment) && longest_fragment_length > 0
            then "-m " + longest_fragment_length
            else ""
        }

        LANG=en; sort -k2,2V -k3,3n -k5,5V -k6,6n -o ~{outfile_name} ~{outfile_name}
    >>>

    output {
        File valid_pairs = outfile_name
        File rs_stats = basename(outfile_name, ".validPairs") + ".RSstat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task merge_valid_interactions {
    meta {
        description: "Merge valid pairs files"
        outputs: {
            all_valid_pairs: "Merged valid pairs file",
            on_target: "On target valid pairs file",
            split_interactions: "Split interactions file",
            all_valid_pairs_stats: "Merged valid pairs statistics file",
        }
    }

    parameter_meta {
        interactions: "Valid pairs files to merge"
        prefix: "Prefix for the output files"
        remove_duplicates: "Remove duplicates"
        report_capture: "Report capture"
        allele_specific_snp: "Use allele specific SNP"
        capture_target: "Capture target file"
    }

    input {
        Array[File] interactions
        String prefix
        Boolean remove_duplicates = false
        Boolean report_capture = false
        Boolean allele_specific_snp = false
        String? capture_target

    }

    command <<<
        set -euo pipefail

        if ~{remove_duplicates}
        then
            LANG=en; sort \
                -k2,2V -k3,3n -k5,5V -k6,6n \
                -m ~{sep(" ", interactions)} \
                | awk -F "\t" \
                'BEGIN{c1=0;c2=0;s1=0;s2=0}(c1!=$2 || c2!=$5 || s1!=$3 || s2!=$6){print;c1=$2;c2=$5;s1=$3;s2=$6}' \
                > ~{prefix}.allValidPairs
        else
            cat ~{sep(" ", interactions)} > ~{prefix}.allValidPairs
        fi

        allcount=$(cat ~{sep(" ", interactions)} | wc -l)
        allcount_rmdup=$(cat ~{prefix}.allValidPairs | wc -l)
        ndbup=$(( $allcount - $allcount_rmdup ))

        ## merge stat file
        echo -e "valid_interaction\t"$allcount > ~{prefix}_allValidPairs.mergestat
        echo -e "valid_interaction_rmdup\t"$allcount_rmdup >> ~{prefix}_allValidPairs.mergestat
        awk 'BEGIN{cis=0;trans=0;sr=0;lr=0} $2 == $5{cis=cis+1; d=$6>$3?$6-$3:$3-$6; if (d<=20000){sr=sr+1}else{lr=lr+1}} $2!=$5{trans=trans+1}END{print "trans_interaction\t"trans"\ncis_interaction\t"cis"\ncis_shortRange\t"sr"\ncis_longRange\t"lr}' ~{prefix}.allValidPairs >> ~{prefix}_allValidPairs.mergestat

        if ~{if defined(capture_target) then "true" else "false"}
        then
            python /HiC-Pro_3.0.0/scripts/onTarget.py \
                -i ~{prefix}.allValidPairs \
                -t ~{capture_target} \
                ~{if (report_capture) then "--cis" else ""} \
                -s ~{prefix}_allValidPairs.mergestat \
                -v \
                > ~{prefix}_ontarget.allValidPairs
        fi

        if ~{allele_specific_snp}
        then
            python /HiC-Pro_3.0.0/scripts/split_valid_interactions.py \
                ~{if defined(capture_target)
                    then "-i  ~{prefix}_ontarget.allValidPairs"
                    else "-i ~{prefix}.allValidPairs"
                } \
                -s ~{prefix}_allValidPairs_assplit.stat \
                -v    
        fi
    >>>

    output {
        File all_valid_pairs = prefix + ".allValidPairs"
        File? on_target = prefix + "_ontarget.allValidPairs"
        File? split_interactions = prefix + "_allValidPairs_assplit.stat"
        File all_valid_pairs_stats = prefix + "_allValidPairs.mergestat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task merge_stats {
    meta {
        description: "Merge statistics files"
        outputs: {
            read1_mapping_stats_merged: "Mapping statistics for read 1",
            read2_mapping_stats_merged: "Mapping statistics for read 2",
            valid_pairs_stats_merged: "Valid pairs statistics",
            rs_stats_merged: "RS statistics",
        }
    }

    parameter_meta {
        read1_mapping_stats: "Mapping statistics for read 1"
        read2_mapping_stats: "Mapping statistics for read 2"
        valid_pairs_stats: "Valid pairs statistics"
        rs_stats: "RS statistics"
        prefix: "Prefix for the output files"
    }

    input {
        Array[File] read1_mapping_stats
        Array[File] read2_mapping_stats
        Array[File] valid_pairs_stats
        Array[File] rs_stats
        String prefix
    }

    command <<<
        set -euo pipefail

        # merge read1 mapping stats
        mkdir read1
        ln -s ~{sep(" ", read1_mapping_stats)} read1/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d read1 \
            -p "*.R1*.mapstat" \
            -v \
            > ~{prefix}_R1.mmapstat
        # merge read2 mapping stats
        mkdir read2
        ln -s ~{sep(" ", read2_mapping_stats)} read2/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d read2 \
            -p "*.R2*.mapstat" \
            -v \
            > ~{prefix}_R2.mmapstat
        # merge pairing stats
        mkdir pairs
        ln -s ~{sep(" ", valid_pairs_stats)} pairs/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d pairs \
            -p "*.pairstat" \
            -v \
            > ~{prefix}.mpairstat
        # merge RS stat
        mkdir rsstat
        ln -s ~{sep(" ", rs_stats)} rsstat/
        python /HiC-Pro_3.0.0/scripts/merge_statfiles.py \
            -d rsstat \
            -p "*.RSstat" \
            -v \
            > ~{prefix}.mRSstat
    >>>

    output {
        File read1_mapping_stats_merged = prefix + "_R1.mmapstat"
        File read2_mapping_stats_merged = prefix + "_R2.mmapstat"
        File valid_pairs_stats_merged = prefix + ".mpairstat"
        File rs_stats_merged = prefix + ".mRSstat"
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
    }
}

task build_raw_maps {
    meta {
        description: "Build raw maps from valid pairs for all specified resolutions"
        outputs: {
            contact_counts: "Contact counts files for each resolution"
        }
    }

    parameter_meta {
        chromsizes: "Chromosome sizes file"
        hic_file: "Hi-C file"
        bin_sizes: "Bin sizes for each contact counts file"
        prefix: "Prefix for the output files"
        all_valid: "Use all valid pairs"
        allele_specific_snp: "Use allele specific SNP"
        capture_target: "Capture target file"
        genome_fragment: "Genome fragment file"
        matrix_format: "Matrix format"
    }

    input {
        File chromsizes
        File hic_file
        Array[Int] bin_sizes
        String prefix
        Boolean allele_specific_snp = false
        Boolean capture_target = false
        File? genome_fragment
        String matrix_format = "upper"
    }

    command <<<
        set -euo pipefail

        for bin_size in ~{sep(" ", bin_sizes)}
        do
            cat ~{hic_file} \
            | /HiC-Pro_3.0.0/scripts/build_matrix \
                --matrix-format ~{matrix_format} \
                ~{if (length(bin_sizes) > 1)
                    then "--binsize $bin_size"
                    else "--binfile ~{select_first([genome_fragment, ""])}"
                } \
                --chrsizes ~{chromsizes} \
                --ifile /dev/stdin \
                --oprefix ~{prefix}_${bin_size}

        done
    >>>

    output {
        Array[File] contact_counts = glob(prefix + "_*.matrix")
        Array[File] contact_counts_bed = glob(prefix + "_*.bed")
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
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

task ice_normalization {
    meta {
        description: "Normalize Hi-C data using the ICE algorithm to correct several sources of bias"
        external_help: "https://nservant.github.io/HiC-Pro/MANUAL.html"
    }

    parameter_meta {
        contact_counts: "Contact counts files to normalize"
        bin_sizes: "Bin sizes for each contact counts file"
        prefix: "Prefix for the output files"
        filtering_percentage: "Percentage of contacts to filter out"
        remove_all_zeroes_loci: "Remove all zero loci"
        filter_low_counts_percentage: "Filter low counts percentage"
        filter_high_counts_percentage: "Filter high counts percentage"
        precision: "Precision of the normalization"
        max_iterations: "Maximum number of iterations"
    }

    input {
        Array[File] contact_counts
        Array[Int] bin_sizes
        String prefix
        Float? filtering_percentage
        Boolean remove_all_zeroes_loci = false
        Float filter_low_counts_percentage = 0.02
        Float filter_high_counts_percentage = 0.0
        Float precision = 0.1
        Int max_iterations = 100
    }

    command <<<
        for map in ~{sep(" ", contact_counts)}
        do
            name=$(basename $map)
            for bin in ~{sep(" ", bin_sizes)}
            do
                ice \
                    --results_filename ${name}_${bin}_iced.matrix \
                    --filter_low_counts_perc ~{filter_low_counts_percentage} \
                    --filter_high_counts_perc ~{filter_high_counts_percentage} \
                    --max_iter ~{max_iterations} \
                    --eps ~{precision} \
                    --remove-all-zeros-loci \
                    --output-bias 1 \
                    $map
            done
        done
    >>>

    output {
        Array[File] iced_matrices = glob(prefix + "*_iced.matrix")
        Array[File] biases = glob(prefix + "*_iced.matrix.biases")
    }

    runtime {
        container: "nservant/hicpro:3.0.0"
        maxRetries: 1
        cpu: 1
        memory: "4 GB"
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

task qcreport {
    meta {
        description: "Write QC report for Hi-C experiment"
        outputs: {
            qc_report: "QC report file"
        }
    }

    parameter_meta {

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