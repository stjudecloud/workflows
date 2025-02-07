version 1.1

task extract_promoters {
    meta {
        description: "Extract promoters from a GTF file"
        help: "Find all transcripts, excluding the mitochondrial chromosome. Add a 2000bp buffer to the transcription start site to include promoter. Extract gene_id and transcript_id ($10 and $12). Subtract the suffix from the transcript ID. Convert to BED format"
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
        container: "ghcr.io/stjudecloud/util:2.0.0"
    }
}

task extract_genes {
    meta {
        description: "Extract genes from a GTF file"
        help: "Find all transcripts, excluding the mitochondrial chromosome. Extract gene_id and transcript_id ($10 and $12). Subtract the suffix from the transcript ID. Convert to BED format"
        outputs: {
            genes: "Gene regions in BED format"
        }
    }

    parameter_meta {
        annotation: "GTF (optionally gzip compressed) file containing gene annotations"
        prefix: "Prefix for output regions file. The extension `.bed` will be added."
    }

    input {
        File annotation
        String prefix = "gene_regions"
    }

    String base = basename(annotation, ".gz")

    #@ except: LineWidth
    command <<<
        set -euo pipefail

        gunzip -c ~{annotation} > ~{base} \
           || ln -sf ~{annotation} ~{base}

        sed -i 's/ /\t/g' ~{base}

        awk -F\\t '{ if($3 == "transcript" && $1 !~ "chrM") 
            print $1 "++" $4 "++" $5 "++" substr($12,2,length($12)-5) "|" substr($10,2,length($10)-3)}' \
            ~{base} \
            | grep -v "unassigned" \
            | sort \
            | uniq \
            | sed s/\+\+/\\t/g \
            | sort -k1,1 -k2,2n > ~{prefix}.bed
    >>>

    output {
        File genes = prefix + ".bed"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        maxRetries: 1
        container: "ghcr.io/stjudecloud/util:2.0.0"
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
        container: "ghcr.io/stjudecloud/hilow:branch-hic_workflow-1.0.0"
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

        python /usr/local/bin/filter_hic.py \
            --all_valid_pairs ~{all_valid_pairs} \
            --filter_pair filter.pair \
            --prefix ~{prefix}

        all=$(wc -l ~{all_valid_pairs} |cut -d " " -f 1)
        filtered=$(wc -l ~{prefix}.allValidPairs.removed | cut -d " " -f 1)
        percentage=$(echo "$filtered/$all*100" | bc -l)
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
        container: "ghcr.io/stjudecloud/hilow:branch-hic_workflow-1.0.0"
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
        container: "ghcr.io/stjudecloud/juicertools:branch-hic_workflow-1.6.2"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        memory: "18 GB"
        maxRetries: 1
    }
}
