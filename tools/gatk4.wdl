## [Homepage](https://software.broadinstitute.org/gatk)
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

task split_n_cigar_reads {
    meta {
        description: "Splits reads that contain Ns in their CIGAR strings into multiple reads."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads"
        outputs: {
            split_n_reads_bam: "BAM file with reads split at N CIGAR elements and updated CIGAR strings.",
            split_n_reads_bam_index: "Index file for the split BAM",
            split_n_reads_bam_md5: "MD5 checksum for the split BAM"
        }
    }

    parameter_meta  {
        bam: "Input BAM format file to with unsplit reads containing Ns in their CIGAR strings."
        bam_index: "BAM index file corresponding to the input BAM"
        fasta: "Reference genome in FASTA format. Must be uncompressed."
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        interval_list: "Interval list indicating regions in which to split reads"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File bam
        File bam_index
        File interval_list
        File fasta
        File fasta_index
        File dict
        String prefix = basename(bam, ".bam") + ".split"
        Int modify_disk_size_gb = 0
        Int ncpu = 8
    }

    Int disk_size_gb = ceil(size(bam, "GB") + 1) * 3 + ceil(size(fasta, "GB")) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        gatk \
            SplitNCigarReads \
            -R ~{fasta} \
            -I ~{bam} \
            -O ~{prefix}.bam \
            -OBM true
       # GATK is unreasonable and uses the plain ".bai" suffix.
       mv ~{prefix}.bai ~{prefix}.bam.bai
    >>>

    output {
        File split_n_reads_bam = "~{prefix}.bam"
        File split_n_reads_bam_index = "~{prefix}.bam.bai"
        File split_n_reads_bam_md5 = "~{prefix}.bam.md5"
    }

    runtime {
        cpu: ncpu
        memory: "24 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task base_recalibrator {
    meta {
        description: "Generates recalibration report for base quality score recalibration."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036897372-BaseRecalibratorSpark-BETA"
        outputs: {
            recalibration_report: "Recalibration report file"
        }
    }

    parameter_meta  {
        bam: "Input BAM format file on which to recabilbrate base quality scores"
        bam_index: "BAM index file corresponding to the input BAM"
        outfile_name: "Name for the output recalibration report."
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        known_indels_sites_VCFs: "List of VCF files containing known indels"
        known_indels_sites_indices: "List of VCF index files corresponding to the VCF files in `known_indels_sites_VCFs`"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
    }

    input {
        File bam
        File bam_index
        String outfile_name = basename(bam, ".bam") + ".recal.txt"
        File fasta
        File fasta_index
        File dict
        File dbSNP_vcf
        File dbSNP_vcf_index
        Array[File] known_indels_sites_VCFs
        Array[File] known_indels_sites_indices
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 4
        Boolean use_original_quality_scores = true
    }

    Int disk_size_gb = ceil(size(bam, "GB") + 1) * 3 + ceil(size(fasta, "GB")) + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        gatk \
            --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4000m -Xmx~{java_heap_size}g" \
            BaseRecalibratorSpark \
            -R ~{fasta} \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{outfile_name} \
            -known-sites ~{dbSNP_vcf} \
            -known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
            --spark-master local[~{ncpu}]
    >>>

    output {
        File recalibration_report = "~{outfile_name}"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task apply_bqsr {
    meta {
        description: "Applies base quality score recalibration to a BAM file."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360040097972-ApplyBQSRSpark-BETA"
        outputs: {
            recalibrated_bam: "Recalibrated BAM file",
            recalibrated_bam_index: "Index file for the recalibrated BAM"
        }

    }

    parameter_meta  {
        bam: "Input BAM format file on which to apply base quality score recalibration"
        bam_index: "BAM index file corresponding to the input BAM"
        recalibration_report: "Recalibration report file"
        prefix: "Prefix for the output recalibrated BAM. The extension `.bqsr.bam` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
    }

    input {
        File bam
        File bam_index
        File recalibration_report
        String prefix = basename(bam, ".bam")
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 4
        Boolean use_original_quality_scores = false
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 2) + 30 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        gatk \
            --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx~{java_heap_size}g" \
            ApplyBQSRSpark \
            --spark-master local[~{ncpu}] \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{prefix}.bqsr.bam \
            --bqsr-recal-file ~{recalibration_report}
    >>>

    output {
        File recalibrated_bam = "~{prefix}.bqsr.bam"
        File recalibrated_bam_index = "~{prefix}.bqsr.bam.bai"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task haplotype_caller {
    meta {
        description: "Calls germline SNPs and indels via local re-assembly of haplotypes."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller"
        outputs: {
            vcf: "VCF file containing called variants",
            vcf_index: "Index file for the VCF"
        }
    }

    parameter_meta  {
        bam: "Input BAM format file on which to call variants"
        bam_index: "BAM index file corresponding to the input BAM"
        interval_list: {
            description: "Interval list indicating regions in which to call variants",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists"
        }
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        prefix: "Prefix for the output VCF. The extension `.vcf.gz` will be added."
        stand_call_conf: {
            description: "Minimum confidence threshold for calling variants",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller#--standard-min-confidence-threshold-for-calling"
        }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
        use_soft_clipped_bases: "Use soft clipped bases in variant calling. Default is to ignore soft clipped bases."
    }

    input {
        File bam
        File bam_index
        File interval_list
        File fasta
        File fasta_index
        File dict
        File dbSNP_vcf
        File dbSNP_vcf_index
        String prefix = basename(bam, ".bam")
        Int stand_call_conf = 20
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 4
        Boolean use_soft_clipped_bases = false
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 2) + 30 + ceil(size(fasta, "GB")) + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        gatk \
           --java-options "-Xms6000m -Xmx~{java_heap_size}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R ~{fasta} \
            -I ~{bam} \
            -L ~{interval_list} \
            -O ~{prefix}.vcf.gz \
            ~{if use_soft_clipped_bases then "" else "--dont-use-soft-clipped-bases"} \
            --standard-min-confidence-threshold-for-calling ~{stand_call_conf} \
            --dbsnp ~{dbSNP_vcf}
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task variant_filtration {
    meta {
        description: "Filters variants based on specified criteria."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration"
        outputs: {
            vcf_filtered: "Filtered VCF file",
            vcf_filtered_index: "Index file for the filtered VCF"
        }
    }

    parameter_meta  {
        vcf: "Input VCF format file to filter"
        vcf_index: "VCF index file corresponding to the input VCF"
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        prefix: "Prefix for the output filtered VCF. The extension `.filtered.vcf.gz` will be added."
        filter_names: {
            description: "Names of the filters to apply",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration#--filter-name"
        }
        filter_expressions: {
            description: "Expressions for the filters",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration#--filter-expression"
        }
        cluster: "Number of SNPs that must be present in a window to filter"
        window: "Size of the window (in bases) for filtering"
        modify_disk_size_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File vcf
        File vcf_index
        File fasta
        File fasta_index
        File dict
        String prefix = basename(vcf, ".vcf.gz")
        Array[String] filter_names = ["FS", "QD"]
        Array[String] filter_expressions = ["FS > 30.0", "QD < 2.0"]
        Int cluster = 3
        Int window = 35
        Int modify_disk_size_gb = 0
        Int ncpu = 1
    }

    Int disk_size_gb = ceil(size(vcf, "GB") * 2) + 30 + modify_disk_size_gb

    command <<<
        gatk VariantFiltration \
            --R ~{fasta} \
                --V ~{vcf} \
                --window ~{window} \
                --cluster ~{cluster} \
                 ~{sep(' ', prefix('--filter-name ', filter_names))} \
                 ~{sep(' ', prefix('--filter-expression ', squote(filter_expressions)))} \
                -O ~{prefix}.filtered.vcf.gz
    >>>

    output {
        File vcf_filtered = "~{prefix}.filtered.vcf.gz"
        File vcf_filtered_index = "~{prefix}.filtered.vcf.gz.tbi"
    }

    runtime {
        cpu: ncpu
        memory: "15 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}
