## [Homepage](https://software.broadinstitute.org/gatk)

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
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
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
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 8
    }

    Int disk_size_gb = ceil(size(bam, "GB") + 1) * 3 + ceil(size(fasta, "GB")) + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        gatk \
            --java-options "-Xms4000m -Xmx~{java_heap_size}g" \
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
        memory: "~{memory_gb} GB"
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
        known_indels_sites_vcfs: "List of VCF files containing known indels"
        known_indels_sites_indices: "List of VCF index files corresponding to the VCF files in `known_indels_sites_vcfs`"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
    }

    input {
        File bam
        File bam_index
        File fasta
        File fasta_index
        File dict
        #@ except: SnakeCase
        File dbSNP_vcf
        #@ except: SnakeCase
        File dbSNP_vcf_index
        Array[File] known_indels_sites_vcfs
        Array[File] known_indels_sites_indices
        String outfile_name = basename(bam, ".bam") + ".recal.txt"
        Boolean use_original_quality_scores = false
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 4
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
            --known-sites ~{dbSNP_vcf} \
            ~{sep(" ", prefix("--known-sites ", known_indels_sites_vcfs))} \
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
        Boolean use_original_quality_scores = false
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 4
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
        #@ except: SnakeCase
        File dbSNP_vcf
        #@ except: SnakeCase
        File dbSNP_vcf_index
        String prefix = basename(bam, ".bam")
        Boolean use_soft_clipped_bases = false
        Int stand_call_conf = 20
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 4
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
        Array[String] filter_names = ["FS", "QD"]
        Array[String] filter_expressions = ["FS > 30.0", "QD < 2.0"]
        String prefix = basename(vcf, ".vcf.gz")
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
                 ~{sep(" ", prefix("--filter-name ", filter_names))} \
                 ~{sep(" ", prefix("--filter-expression ", squote(filter_expressions)))} \
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

task mark_duplicates_spark {
     meta {
        description: "Marks duplicate reads in the input BAM file using GATK's Spark implementation of Picard's MarkDuplicates."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/13832682540699-MarkDuplicatesSpark"
        outputs: {
            duplicate_marked_bam: "The input BAM with computationally determined duplicates marked.",
            duplicate_marked_bam_index: "The `.bai` BAM index file associated with `duplicate_marked_bam`",
            mark_duplicates_metrics: {
                description: "The METRICS_FILE result of `picard MarkDuplicates`",
                external_help: "http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics",
            }
        }
    }

    parameter_meta {
        bam: "Input BAM format file in which to mark duplicates"
        prefix: "Prefix for the MarkDuplicates result files. The extensions `.bam`, `.bam.bai`, and `.metrics.txt` will be added."
        duplicate_scoring_strategy: {
            description: "Strategy for scoring duplicates.",
            choices: [
                "SUM_OF_BASE_QUALITIES",
                "TOTAL_MAPPED_REFERENCE_LENGTH",
                "RANDOM"
            ]
        }
        read_name_regex: "Regular expression for extracting tile names, x coordinates, and y coordinates from read names. The default works for typical Illumina read names."
        tagging_policy: {
            description: "Tagging policy for the output BAM.",
            choices: [
                "DontTag",
                "OpticalOnly",
                "All"
            ],
        }
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT"
            ],
            tool_default: "STRICT",
        }
        create_bam: {
            description: "Enable BAM creation (true)? Or only output MarkDuplicates metrics (false)?",
            common: true
        }
        optical_distance: "Maximum distance between read coordinates to consider them optical duplicates. If `0`, then optical duplicate marking is disabled. Suggested settings of 100 for unpatterned versions of the Illumina platform (e.g. HiSeq) or 2500 for patterned flowcell models (e.g. NovaSeq). Calculation of distance depends on coordinate data embedded in the read names, typically produced by the Illumina sequencing machines. Optical duplicate detection will not work on non-standard names without modifying `read_name_regex`."
        modify_memory_gb: "Add to or subtract from the default memory allocation. Default memory allocation is determined by the size of the input BAM. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".MarkDuplicates"
        String duplicate_scoring_strategy = "SUM_OF_BASE_QUALITIES"
        String read_name_regex = "^[!-9;-?A-~:]+:([!-9;-?A-~]+):([0-9]+):([0-9]+)$"
        String tagging_policy = "All"
        String validation_stringency = "SILENT"
        Boolean create_bam = true
        Int optical_distance = 0
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
        Int ncpu = 4
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = min(ceil(bam_size + 15), 50) + modify_memory_gb
    Int disk_size_gb = (
        (
            if create_bam
            then ceil((bam_size * 2) + 10)
            else ceil(bam_size + 10)
        ) + modify_disk_size_gb
    )

    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        gatk MarkDuplicatesSpark \
            --java-options "-Xmx~{java_heap_size}g" \
            -I ~{bam} \
            -M ~{prefix}.metrics.txt \
            -O ~{if create_bam then prefix + ".bam" else "/dev/null"} \
            --create-output-bam-index ~{create_bam} \
            --read-validation-stringency ~{validation_stringency} \
            --duplicate-scoring-strategy ~{duplicate_scoring_strategy} \
            --read-name-regex '~{
                if (optical_distance > 0) then read_name_regex else "null"
            }' \
            --duplicate-tagging-policy ~{tagging_policy} \
            --optical-duplicate-pixel-distance ~{optical_distance} \
            --spark-master local[~{ncpu}]
    >>>

    output {
        File? duplicate_marked_bam = "~{prefix}.bam"
        File? duplicate_marked_bam_index = "~{prefix}.bam.bai"
        File mark_duplicates_metrics = "~{prefix}.metrics.txt"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"
        maxRetries: 1
    }
}
