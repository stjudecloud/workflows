## # Utilities
##
## This WDL file includes custom scripts to parse and reformat 
## task output as part of a workflow as well as generic tasks.

version 1.0

task download {
    input {
        String url
        String outfilename
        String? md5sum
        Int disk_size_gb = 10
        Int max_retries = 1
    }

    runtime {
        memory: "4 GB"
        disk: disk_size_gb + " GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }

    command <<<
        set -euo pipefail

        wget ~{url} -O ~{outfilename}

        if [ ! -z "~{md5sum}" ]; then
            echo "~{md5sum}  ~{outfilename}" > ~{outfilename}.md5
            md5sum -c ~{outfilename}.md5
        fi
    >>>

    output {
        File outfile = outfilename
    }

    meta {
        author: "Clay McLeod"
        email: "clay.mcleod@stjude.org"
        description: "This WDL task uses wget to download a file from a remote URL to the local filesystem" 
    }

    parameter_meta {
        url: "URL of the file to download"
        outfilename: "Name to use for the output file"
    }
}

task get_read_groups {
    input {
        File bam
        Int max_retries = 1
        Boolean format_for_star = true
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        set -euo pipefail

        if [ "~{format_for_star}" == "true" ]
        then        
            samtools view -H ~{bam} | grep "@RG" \
                | cut -f 2- \
                | sed -e 's/\t/ /g' \
                | awk '{print}' ORS=' , ' \
                | sed 's/ , $//' >> read_groups.txt
        else
            samtools view -H ~{bam} | grep "@RG" >> read_groups.txt
        fi
    >>>

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }

    output { 
        File out = "read_groups.txt"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL task is a utility to get read group information from a BAM file and write it out to as a string" 
    }

    parameter_meta {
        bam: "Input BAM format file to get read groups from"
    }
}

task split_string {
    input {
        String input_string
        String delimiter = " , "
        Int max_retries = 1
        Int disk_size = 1
    }
    command <<<
        set -euo pipefail

        echo ~{input_string} | sed 's/~{delimiter}/\n/g' > output.txt
    >>>

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }

    output {
        File output_file = "output.txt"
        Array[String] out = read_lines("output.txt")
    }
}

task calc_gene_lengths {
    input {
        File gtf
        String outfile = basename(gtf, ".gtf.gz") + ".genelengths.txt"
        Int max_retries = 1
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(gtf_size * 2 + 10)

    command <<<
        set -euo pipefail

        GTF="~{gtf}" OUTFILE="~{outfile}" python - <<END
import os  # lint-check: ignore
import gtfparse  # lint-check: ignore
import numpy as np  # lint-check: ignore

gtf_name = os.environ['GTF']
outfile = open(os.environ['OUTFILE'], 'w')

gtf = gtfparse.read_gtf(gtf_name)

only_genes = gtf[gtf['feature'] == 'gene']
only_exons = gtf[gtf['feature'] == 'exon']
gene_start_offset = {}
gene_end_offset = {}
gene_exon_intersection = {}
gene_total_exon_size = {}
gene_length = {}

for (index, value) in only_genes.iterrows():
    gene_name = value['gene_name']
    start = value['start']
    end = value['end']
    size = end - start
    
    if size <= 0:
        raise RuntimeError("Size of gene is negative!")
        
    gene_start_offset[gene_name] = start
    gene_end_offset[gene_name] = end
    gene_exon_intersection[gene_name] = np.zeros(size)
    gene_total_exon_size[gene_name] = 0
    gene_length[gene_name] = end - start
    
for (index, value) in only_exons.iterrows():
    gene_name = value['gene_name']
    offset = gene_start_offset[gene_name]
    start = value['start'] - offset
    end = value['end'] - offset
    exon_length = end - start
    gene_exon_intersection[gene_name][start:end] = 1
    gene_total_exon_size[gene_name] += exon_length

results = []
print("Gene name\tlength", file=outfile)
for (gene, exonic_intersection) in sorted(gene_exon_intersection.items()):
    length = np.sum(exonic_intersection).astype(int)
    print(f"{gene}\t{length}", file=outfile)
END
    >>>

    runtime {
        memory: "8 GB"
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0'
        maxRetries: max_retries
    }

    output {
        File out = "~{outfile}"
    }
}

task qc_summary {
    input {
        File multiqc_tar_gz
        String outfile = basename(multiqc_tar_gz, ".multiqc.tar.gz") + ".qc_summary.json"
        Int disk_size = 1
        Int max_retries = 1
    }

    String sample_name = basename(multiqc_tar_gz, ".multiqc.tar.gz")

    command <<<
        set -euo pipefail

        tar -xzf "~{multiqc_tar_gz}"
        gen_stats_file=~{sample_name}.multiqc/multiqc_data/multiqc_general_stats.txt

        TOTAL_READS=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-total_reads $gen_stats_file | tail -n 1 | awk '{ printf("%.0f", $1) }')
        PERCENT_ALIGNED=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-percentage_aligned $gen_stats_file | tail -n 1 | awk '{ printf("%.3f", $1) }')
        MEAN_COVERAGE=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-mean_coverage $gen_stats_file | tail -n 1 | awk '{ printf("%.3f", $1) }')
        INSERT_SIZE=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-median_insert_size $gen_stats_file | tail -n 1)
        MEAN_GC_CONTENT=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-avg_gc $gen_stats_file | tail -n 1 | awk '{ printf("%.3f", $1) }')
        READ_LENGTH=$(csvcut -t -c ngsderive_mqc-generalstats-ngsderive-consensusreadlength $gen_stats_file | tail -n 1)
        PLATFORM=$(csvcut -t -c ngsderive_mqc-generalstats-ngsderive-instrument $gen_stats_file | tail -n 1)

        THIRTYX_PERCENT=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-30_x_pc $gen_stats_file | tail -n 1)
        DUP_PERCENT=$(csvcut -t -c FastQC_mqc-generalstats-fastqc-percent_duplicates $gen_stats_file | tail -n 1)

        STRANDEDNESS=$({ csvcut -t -c ngsderive_mqc-generalstats-ngsderive-predicted $gen_stats_file || echo "Not Applicable" ; } | tail -n 1)

        jq -n \
            --arg TOTAL_READS "$TOTAL_READS" \
            --arg PERCENT_ALIGNED "$PERCENT_ALIGNED" \
            --arg MEAN_COVERAGE "$MEAN_COVERAGE" \
            --arg INSERT_SIZE "$INSERT_SIZE" \
            --arg MEAN_GC_CONTENT "$MEAN_GC_CONTENT" \
            --arg READ_LENGTH "$READ_LENGTH" \
            --arg PLATFORM "$PLATFORM" \
            --arg THIRTYX_PERCENT "$THIRTYX_PERCENT" \
            --arg DUP_PERCENT "$DUP_PERCENT" \
            --arg STRANDEDNESS "$STRANDEDNESS" \
            '{
                total_reads: ($TOTAL_READS | tonumber),
                percent_aligned: ($PERCENT_ALIGNED | tonumber),
                mean_coverage: ($MEAN_COVERAGE | tonumber),
                median_insert_size: ($INSERT_SIZE | tonumber),
                mean_percent_GC: ($MEAN_GC_CONTENT | tonumber),
                inferred_read_length: ($READ_LENGTH | tonumber),
                inferred_sequencing_platform: $PLATFORM,
                percent_thirtyX_coverage: ($THIRTYX_PERCENT | tonumber),
                percent_duplicate: ($DUP_PERCENT | tonumber),
                inferred_strandedness: $STRANDEDNESS 
            }' > ~{outfile}
    >>>

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }

    output {
        File out = "~{outfile}"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task pulls out keys metrics that can provide a high level overview of the sample, without needing to examine the entire MultiQC report. Currently, these key metrics come from Qualimap and ngsderive." 
    }
}

task compression_integrity {
    input {
        File bam
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size + 10)

    command {
        bgzip -t ~{bam}
    }

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task add_to_bam_header {
    input {
        File input_bam
        String additional_header
        String output_bam_name = basename(input_bam, ".bam") + ".reheader.bam"
        Int max_retries = 1
    }

    Float bam_size = size(input_bam, "GiB")
    Int disk_size = ceil(bam_size + 1)

    command <<<
        samtools view -H ~{input_bam} > header.sam
        echo "~{additional_header}" >> header.sam
        samtools reheader -P header.sam ~{input_bam} > ~{output_bam_name}
    >>>

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }

    output {
        File output_file = "header.sam"
        Array[String] out = read_lines("header.sam")
        File output_bam = output_bam_name
    }
}

task unpack_tarball {
    input {
        File tarball
        Int extra_disk_gb = 0
        Int max_retries = 1
    }

    Float tarball_size = size(tarball, "GiB")
    Int disk_size = ceil(tarball_size * 8) + extra_disk_gb

    command <<<
        set -euo pipefail

        mkdir unpacked_tarball
        tar -C unpacked_tarball -xzf ~{tarball}
        find unpacked_tarball/ -type f > file_list.txt
    >>>

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/util:1.2.0'
        maxRetries: max_retries
    }

    output {
        Array[File] tarball_contents = read_lines("file_list.txt")
    }
}