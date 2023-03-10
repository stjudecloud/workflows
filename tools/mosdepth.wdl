## # Mosdepth
##
## This WDL tool wraps the [mosdepth tool](https://github.com/brentp/mosdepth).

version 1.0

task make_coverage_regions_beds {
    input {
        File gtf
        Int max_retries = 1
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(gtf_size * 3)

    command <<<
        set -euo pipefail

        BED=$(basename ~{gtf} '.gz').bed
        gunzip -c ~{gtf} | gtf2bed > "$BED"

        EXON=$(basename ~{gtf} '.gz').exon.bed
        awk '/\texon\t/ {print $1 "\t" $2 "\t" $3}' "$BED" > "$EXON"

        CDS=$(basename ~{gtf} '.gz').CDS.bed
        awk '/\tCDS\t/ {print $1 "\t" $2 "\t" $3}' "$BED" > "$CDS"
    >>>

    runtime {
        disk: disk_size + " GB"
        memory: "4 GB"
        docker: 'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0'
        maxRetries: max_retries
    }

    output {
        File exon_bed = basename(gtf, '.gz') + ".exon.bed"
        File CDS_bed = basename(gtf, '.gz') + ".CDS.bed"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL task takes in a GTF file, converts it to BED, then filters it down to two 3 column BED files: one of only 'exons', one of only 'CDS' regions"
    }
}

task coverage {
    input {
        File bam
        File bai
        File? coverage_bed
        String prefix = basename(bam, '.bam')
        Int min_mapping_quality = 20
        Boolean use_fast_mode = true
        Int memory_gb = 8 
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil(bam_size + 5)

    command {
        set -euo pipefail

        mv ~{bai} ~{bam}.bai || true
        mosdepth \
            -n \
            ~{if defined(coverage_bed) then "-b" else ""} ~{coverage_bed} \
            -Q ~{min_mapping_quality} \
            ~{if (use_fast_mode) then "-x" else ""} \
            ~{prefix} \
            ~{bam}
    }

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'ghcr.io/stjudecloud/mosdepth:1.0.0'
        maxRetries: max_retries
    }

    output {
        File summary = prefix + ".mosdepth.summary.txt"
        File global_dist = prefix + ".mosdepth.global.dist.txt"
        File? region_dist = prefix + ".mosdepth.region.dist.txt"
    }

    meta {
        author: "Andrew Frantz"
        email: "andrew.frantz@stjude.org"
        description: "This WDL tool wraps the Mosdepth tool for calculating coverage"
    }
}
