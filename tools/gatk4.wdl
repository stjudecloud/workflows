## [Homepage](https://software.broadinstitute.org/gatk)
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

task split_n_cigar_reads {
    meta {

    }

    parameter_meta  {
        bam: "Input BAM format file to <brief description of task>"
        bam_index: "BAM index file corresponding to the input BAM"
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        interval_list: "Interval list for splitting reads"
        prefix: "Prefix for the BAM file. The extension `.split.bam` will be added."
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
        gatk \
            SplitNCigarReads \
            -R ~{fasta} \
            -I ~{bam} \
            -O ~{prefix}.bam \
            -OBM true
    >>>

    output {
        File split_bam = "~{prefix}.bam"
        File split_bam_index = "~{prefix}.bam.bai"
        File split_bam_md5 = "~{prefix}.bam.md5"
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task base_recalibrator {
    meta {

    }

    parameter_meta  {
        bam: "Input BAM format file to <brief description of task>"
        bam_index: "BAM index file corresponding to the input BAM"
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        known_indels_sites_VCFs: "List of VCF files containing known indels"
        known_indels_sites_indices: "List of VCF index files corresponding to the VCF files in `known_indels_sites_VCFs`"
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
        prefix: "Prefix for the BAM file. The extension `.split.bam` will be added."
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
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R ~{fasta} \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{prefix}.txt \
            -known-sites ~{dbSNP_vcf} \
            -known-sites ~{sep=" --known-sites " known_indels_sites_VCFs}
    >>>

    output {
        File recalibration_report = "~{prefix}.txt"
    }

    runtime {
        memory: "6 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task apply_bqsr {
    meta {

    }

    parameter_meta  {

    }

    input {
        File bam
        File bam_index
        File recalibration_report
        File dict
        File fasta
        File fasta_index
        String prefix = basename(bam, ".bam") + ".bqsr"
        Int modify_disk_size_gb = 0
        Boolean use_original_quality_scores = true
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 4) + 30 + ceil(size(fasta, "GB")) + modify_disk_size_gb

    command <<<
        gatk \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R ~{fasta} \
            -I ~{bam} \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O ~{prefix}.bam \
            --bqsr-recal-file ~{recalibration_report}
    >>>

    output {
        File recalibrated_bam = "~{prefix}.bam"
        File recalibrated_bam_index = "~{prefix}.bai"
    }

    runtime {
        memory: "3.5 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}



task haplotype_caller {
    meta {

    }

    parameter_meta  {

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
		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
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
        memory: "6.5 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}

task variant_filtration {
    meta {

    }

    parameter_meta  {

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
            ~{sep(' ', prefix('--filter-expression ', filter_expression))} \
			-O ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_filtered = "~{prefix}.vcf.gz"
        File vcf_filtered_index = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        memory: "3 GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}
