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
            split_bam: "BAM file with reads split a N CIGAR elements and updated CIGAR strings."
            split_bam_index: "Index file for the split BAM"
            split_bam_md5: "MD5 checksum for the split BAM"
        }
    }

    parameter_meta  {
        bam: "Input BAM format file to with unsplit reads containg Ns in their CIGAR strings."
        bam_index: "BAM index file corresponding to the input BAM"
        fasta: "Reference genome in FASTA format. Must be uncompressed."
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        interval_list: "Interval list indicating regions in which to split reads"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
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
    }

    Int disk_size_gb = ceil(size(bam, "GB") + 1) * 5 + ceil(size(fasta, "GB")) + modify_disk_size_gb

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
        File split_bam = "~{prefix}.bam"
        File split_bam_index = "~{prefix}.bam.bai"
        File split_bam_md5 = "~{prefix}.bam.md5"
    }

    runtime {
        cpu: 8
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
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        known_indels_sites_VCFs: "List of VCF files containing known indels"
        known_indels_sites_indices: "List of VCF index files corresponding to the VCF files in `known_indels_sites_VCFs`"
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
        prefix: "Prefix for the output recalibration report. The extension `.txt` will be added."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File bam_index
        String prefix = basename(bam, ".bam") + ".recal"
        File dict
        File fasta
        File fasta_index
        File dbSNP_vcf
        File dbSNP_vcf_index
        Array[File] known_indels_sites_VCFs
        Array[File] known_indels_sites_indices
        Int modify_disk_size_gb = 0
        Boolean use_original_quality_scores = true
    }

    Int disk_size_gb = ceil(size(bam, "GB") + 1) * 3 + ceil(size(fasta, "GB")) + modify_disk_size_gb

    command <<<
        gatk \
            --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibratorSpark \
            -R ~{fasta} \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{prefix}.txt \
            -known-sites ~{dbSNP_vcf} \
            -known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
            --spark-master local[4]
    >>>

    output {
        File recalibration_report = "~{prefix}.txt"
    }

    runtime {
        cpu: 4
        memory: "25 GB"
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
            recalibrated_bam: "Recalibrated BAM file"
            recalibrated_bam_index: "Index file for the recalibrated BAM"
        }

    }

    parameter_meta  {
        bam: "Input BAM format file on which to apply base quality score recalibration"
        bam_index: "BAM index file corresponding to the input BAM"
        recalibration_report: "Recalibration report file"
        prefix: "Prefix for the output recalibrated BAM. The extension `.bam` will be added."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
    }

    input {
        File bam
        File bam_index
        File recalibration_report
        String prefix = basename(bam, ".bam") + ".bqsr"
        Int modify_disk_size_gb = 0
        Boolean use_original_quality_scores = true
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 4) + 30 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        gatk \
            --java-options "-XX:+PrintFlagsFinal \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSRSpark \
            --spark-master local[4] \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{prefix}.bam \
            --bqsr-recal-file ~{recalibration_report}
       # GATK is unreasonable and uses the plain ".bai" suffix.
       #mv ~{prefix}.bai ~{prefix}.bam.bai
    >>>

    output {
        File recalibrated_bam = "~{prefix}.bam"
        File recalibrated_bam_index = "~{prefix}.bam.bai"
    }

    runtime {
        cpu: 4
        memory: "25 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task haplotype_caller {
    meta {
        description: "Calls germline SNPs and indels via local re-assembly of haplotypes."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037433931-HaplotypeCallerSpark-BETA"
        outputs: {
            vcf: "VCF file containing called variants"
            vcf_index: "Index file for the VCF"
        }
    }

    parameter_meta  {
        bam: "Input BAM format file on which to call variants"
        bam_index: "BAM index file corresponding to the input BAM"
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        interval_list: "Interval list indicating regions in which to call variants"
        prefix: "Prefix for the output VCF. The extension `.vcf.gz` will be added."
        stand_call_conf: "Minimum confidence threshold for calling variants"
        modify_disk_size_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        use_soft_clipped_bases: "Use soft clipped bases in variant calling. Default is to ignore soft clipped bases."
    }

    input {
        File bam
        File bam_index
        File interval_list
        File dict
        File fasta
        File fasta_index
        File dbSNP_vcf
        File dbSNP_vcf_index
        String prefix = basename(bam, ".bam")
        Int stand_call_conf = 20
        Int modify_disk_size_gb = 0
        Boolean use_soft_clipped_bases = false
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 2) + 30 + ceil(size(fasta, "GB")) + modify_disk_size_gb

    command <<<
		    gatk \
           --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCallerSpark \
            -R ~{fasta} \
            -I ~{bam} \
            -L ~{interval_list} \
            -O ~{prefix}.vcf.gz \
            ~{if use_soft_clipped_bases then "" else "--dont-use-soft-clipped-bases"} \
            --standard-min-confidence-threshold-for-calling ~{stand_call_conf} \
            --dbsnp ~{dbSNP_vcf} \
            --spark-master local[4]
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        cpu: 4
        memory: "25 GB"
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
            vcf_filtered: "Filtered VCF file"
            vcf_filtered_index: "Index file for the filtered VCF"
        }
    }

    parameter_meta  {
        vcf: "Input VCF format file to filter"
        vcf_index: "VCF index file corresponding to the input VCF"
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        prefix: "Prefix for the output filtered VCF. The extension `.vcf.gz` will be added."
        filter_name: "Name of the filter to apply"
        filter_expression: "Expression for the filter"
        cluster: "Number of SNPs that must be present in a window to filter"
        window: "Size of the window for filtering"
        modify_disk_size_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
    }

    input {
        File vcf
        File vcf_index
        File dict
        File fasta
        File fasta_index
        String prefix = basename(vcf, ".vcf.gz") + ".filtered"
        Array[String] filter_name = ["FS", "QD"]
        Array[String] filter_expression = ["FS > 30.0", "QD < 2.0"]
        Int cluster = 3
        Int window = 35
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(vcf, "GB") * 2) + 30 + modify_disk_size_gb

    command <<<
        gatk \
            VariantFiltration \
                --R ~{fasta} \
                --V ~{vcf} \
                --window ~{window} \
                --cluster ~{cluster} \
                 ~{sep(' ', prefix('--filter-name ', filter_name))} \
                 ~{sep(' ', prefix('--filter-expression ', squote(filter_expression)))} \
                -O ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_filtered = "~{prefix}.vcf.gz"
        File vcf_filtered_index = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        cpu: 1
        memory: "15 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}
