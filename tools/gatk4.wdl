## [Homepage](https://software.broadinstitute.org/gatk)

version 1.3

enum ref_confidence {
    NONE,
    GVCF,
    BP_RESOLUTION
}

enum VariantMode {
    SNP,
    INDEL,
    BOTH
}

struct Resource {
    meta {
        description: "Struct representing a resource to be used in GATK VariantRecalibrator. This includes both the VCF file for the resource and metadata about how to use the resource in building the recalibration model."
    }

    parameter_meta {
        name: "Name of the resource. This is an arbitrary string that is used to identify the resource in the recalibration model."
        known: "Boolean indicating whether the resource is a known set of variants."
        training: "Boolean indicating whether the resource is a training set of variants to be used for building the recalibration model."
        truth: "Boolean indicating whether the resource is a truth set of variants to be used for building the recalibration model."
        prior: "Prior probability for the resource."
        vcf: "VCF file for the resource."
        vcf_index: "Index file for the VCF."
    }

    String name
    Boolean known
    Boolean training
    Boolean truth
    Float prior
    File vcf
    File vcf_index
}

task resource_to_string {
    meta {
        description: "Converts a resource struct to a string representation for use in GATK commands."
        outputs: {
            res_string: "String representation of the input resource struct formatted for GATK VariantRecalibrator.",
            res_vcf: "VCF file from the input resource struct",
            res_vcf_index: "Index file for the VCF from the input resource struct",
        }
    }

    parameter_meta {
        res: "Resource struct containing information about a resource to be used in GATK VariantRecalibrator."
    }

    input {
        Resource res
    }

    command <<<
        echo "~{res.name},known=~{res.known},training=~{res.training},truth=~{res.truth},prior=~{res.prior} ~{basename(res.vcf)}"
    >>>

    output {
        String res_string = read_string(stdout())
        File res_vcf = res.vcf
        File res_vcf_index = res.vcf_index
    }

    requirements {
        container: "ghcr.io/stjudecloud/util:3.0.3"
        maxRetries: 1
    }
}

task split_n_cigar_reads {
    meta {
        description: "Splits reads that contain Ns in their CIGAR strings into multiple reads."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads"
        outputs: {
            split_n_reads_bam: "BAM file with reads split at N CIGAR elements and updated CIGAR strings.",
            split_n_reads_bam_index: "Index file for the split BAM",
            split_n_reads_bam_md5: "MD5 checksum for the split BAM",
        }
    }

    parameter_meta  {
        bam: "Input BAM format file to with unsplit reads containing Ns in their CIGAR strings."
        bam_index: "BAM index file corresponding to the input BAM"
        fasta: "Reference genome in FASTA format. Must be uncompressed."
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File bam
        File bam_index
        File fasta
        File fasta_index
        File dict
        String prefix = basename(bam, ".bam") + ".split"
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 8
    }

    Int disk_size_gb = ceil(size(bam, "GB") + 1) * 3
        + ceil(size(fasta, "GB"))
        + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        gatk \
            --java-options "-Xms4000m -Xmx~{java_heap_size}g" \
            SplitNCigarReads \
            -R "~{fasta}" \
            -I "~{bam}" \
            -O "~{prefix}.bam" \
            -OBM true
       # GATK is unreasonable and uses the plain ".bai" suffix.
       mv "~{prefix}.bai" "~{prefix}.bam.bai"
    >>>

    output {
        File split_n_reads_bam = "~{prefix}.bam"
        File split_n_reads_bam_index = "~{prefix}.bam.bai"
        File split_n_reads_bam_md5 = "~{prefix}.bam.md5"
    }

    requirements {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
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
        known_indels_sites_vcfs: "List of VCF files containing known indels"
        known_indels_sites_indices: "List of VCF index files corresponding to the VCF files in `known_indels_sites_vcfs`"
        outfile_name: "Name for the output recalibration report."
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
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
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int ncpu = 4
        }

    Int disk_size_gb = ceil(size(bam, "GB") + 1) * 3
        + ceil(size(fasta, "GB"))
        + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        ref_fasta=~{sub(basename(fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{fasta}" > "$ref_fasta.fa" \
            || ln -sf "~{fasta}" "$ref_fasta.fa"
        ln -sf "~{fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{dict}" "$ref_fasta.dict"

        # shellcheck disable=SC2102
        gatk \
            --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4000m -Xmx~{java_heap_size}g" \
            BaseRecalibratorSpark \
            -R "$ref_fasta.fa" \
            -I "~{bam}" \
            ~{(
                if use_original_quality_scores
                then "--use-original-qualities"
                else ""
            )} \
            -O "~{outfile_name}" \
            --known-sites "~{dbSNP_vcf}" \
            ~{sep(" ", prefix("--known-sites ", squote(known_indels_sites_vcfs)))} \
            --spark-master local[~{ncpu}]

        rm -rf "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict"
    >>>

    output {
        File recalibration_report = "~{outfile_name}"
    }

    requirements {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
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
            recalibrated_bam_index: "Index file for the recalibrated BAM",
        }
    }

    parameter_meta  {
        bam: "Input BAM format file on which to apply base quality score recalibration"
        bam_index: "BAM index file corresponding to the input BAM"
        recalibration_report: "Recalibration report file"
        prefix: "Prefix for the output recalibrated BAM. The extension `.bqsr.bam` will be added."
        use_original_quality_scores: "Use original quality scores from the input BAM. Default is to use recalibrated quality scores."
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File bam
        File bam_index
        File recalibration_report
        String prefix = basename(bam, ".bam")
        Boolean use_original_quality_scores = false
        Int memory_gb = 50
        Int modify_disk_size_gb = 0
        Int ncpu = 4
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 2) + 30 + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    command <<<
        set -euo pipefail

        # shellcheck disable=SC2102
        gatk \
            --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx~{java_heap_size}g" \
            ApplyBQSRSpark \
            --spark-master local[~{ncpu}] \
            -I "~{bam}" \
            ~{if use_original_quality_scores then "--use-original-qualities" else "" } \
            -O "~{prefix}.bqsr.bam" \
            --bqsr-recal-file "~{recalibration_report}"
    >>>

    output {
        File recalibrated_bam = "~{prefix}.bqsr.bam"
        File recalibrated_bam_index = "~{prefix}.bqsr.bam.bai"
    }

    requirements {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
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
            vcf_index: "Index file for the VCF",
        }
    }

    parameter_meta  {
        bam: "Input BAM format file on which to call variants"
        bam_index: "BAM index file corresponding to the input BAM"
        interval_list: {
            description: "Interval list indicating regions in which to call variants",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists",
        }
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        reference_confidence: {
            description: "Reference confidence mode to run HaplotypeCaller in.",
            help: "If `NONE`, HaplotypeCaller will run in default mode and only output variant sites. If `GVCF`, HaplotypeCaller will run in GVCF mode and output both variant and non-variant sites with reference confidence scores. If `BP_RESOLUTION`, HaplotypeCaller will run in GVCF mode but output non-variant sites at base pair resolution instead of block resolution.",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller#--emit-ref-confidence",
        }
        prefix: "Prefix for the output VCF. The extension `.vcf.gz` will be added."
        use_soft_clipped_bases: "Use soft clipped bases in variant calling. Default is to ignore soft clipped bases."
        stand_call_conf: {
            description: "Minimum confidence threshold for calling variants",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller#--standard-min-confidence-threshold-for-calling",
        }
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
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
        ref_confidence reference_confidence = ref_confidence.NONE
        String prefix = basename(bam, ".bam")
        Boolean use_soft_clipped_bases = false
        Int stand_call_conf = 20
        Int memory_gb = 25
        Int modify_disk_size_gb = 0
        Int ncpu = 4
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 2)
        + 30
        + ceil(size(fasta, "GB"))
        + modify_disk_size_gb
    Int java_heap_size = ceil(memory_gb * 0.9)

    String sample = basename(bam)
    String snp_vcf = basename(dbSNP_vcf)
    String snp_vcf_index = basename(dbSNP_vcf_index)

    command <<<
        set -euo pipefail

        ref_fasta=~{sub(basename(fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{fasta}" > "$ref_fasta.fa" \
            || ln -sf "~{fasta}" "$ref_fasta.fa"
        ln -sf "~{fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{dict}" "$ref_fasta.dict"

        ln -sf "~{bam}" "~{sample}"
        ln -sf "~{bam_index}" "~{sample}.bai"

        ln -sf "~{dbSNP_vcf}" "~{snp_vcf}"
        ln -sf "~{dbSNP_vcf_index}" "~{snp_vcf_index}"

        gatk \
           --java-options "-Xms6000m -Xmx~{java_heap_size}g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R "$ref_fasta.fa" \
            -I "~{sample}" \
            -L "~{interval_list}" \
            -O "~{prefix}.vcf.gz" \
            ~{if use_soft_clipped_bases then "" else "--dont-use-soft-clipped-bases"} \
            --standard-min-confidence-threshold-for-calling ~{stand_call_conf} \
            --dbsnp "~{snp_vcf}" \
            --emit-ref-confidence "~{reference_confidence}"

        rm -rf "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict" "~{sample}" "~{sample}.bai" "~{snp_vcf}" "~{snp_vcf_index}"
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    requirements {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        maxRetries: 1
    }
}

task variant_filtration {
    meta {
        description: "Filters variants based on specified criteria."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration"
        outputs: {
            vcf_filtered: "Filtered VCF file",
            vcf_filtered_index: "Index file for the filtered VCF",
        }
    }

    parameter_meta  {
        vcf: "Input VCF format file to filter"
        vcf_index: "VCF index file corresponding to the input VCF"
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        filter_names: {
            description: "Names of the filters to apply",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration#--filter-name",
        }
        filter_expressions: {
            description: "Expressions for the filters",
            external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration#--filter-expression",
        }
        prefix: "Prefix for the output filtered VCF. The extension `.filtered.vcf.gz` will be added."
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
            --R "~{fasta}" \
            --V "~{vcf}" \
            --window ~{window} \
            --cluster ~{cluster} \
            ~{sep(" ", prefix("--filter-name ", quote(filter_names)))} \
            ~{sep(" ", prefix("--filter-expression ", squote(filter_expressions)))} \
            -O "~{prefix}.filtered.vcf.gz"
    >>>

    output {
        File vcf_filtered = "~{prefix}.filtered.vcf.gz"
        File vcf_filtered_index = "~{prefix}.filtered.vcf.gz.tbi"
    }

    requirements {
        cpu: ncpu
        memory: "15 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
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
            },
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
                "RANDOM",
            ],
        }
        read_name_regex: {
            description: "Regular expression for extracting tile names, x coordinates, and y coordinates from read names.",
            help: "The default works for typical Illumina read names.",
        }
        tagging_policy: {
            description: "Tagging policy for the output BAM.",
            choices: [
                "DontTag",
                "OpticalOnly",
                "All",
            ],
        }
        validation_stringency: {
            description: "Validation stringency for parsing the input BAM.",
            choices: [
                "STRICT",
                "LENIENT",
                "SILENT",
            ],
            tool_default: "STRICT",
        }
        create_bam: {
            description: "Enable BAM creation (true)? Or only output MarkDuplicates metrics (false)?",
            group: "Common",
        }
        optical_distance: {
            description:  "Maximum distance between read coordinates to consider them optical duplicates. If `0`, then optical duplicate marking is disabled.",
            help: "Suggested settings of 100 for unpatterned versions of the Illumina platform (e.g. HiSeq) or 2500 for patterned flowcell models (e.g. NovaSeq). Calculation of distance depends on coordinate data embedded in the read names, typically produced by the Illumina sequencing machines.",
            warning: "Optical duplicate detection will not work on non-standard names without modifying `read_name_regex`.",
        }
        modify_memory_gb: "Add to or subtract from the dynamic memory allocation. Default memory allocation is determined by the size of the inputs. Specified in GB."
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

    Float bam_size = size(bam, "GB")
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

        # shellcheck disable=SC2102
        gatk MarkDuplicatesSpark \
            --java-options "-Xmx~{java_heap_size}g" \
            -I "~{bam}" \
            -M "~{prefix}.metrics.txt" \
            -O "~{if create_bam then prefix + ".bam" else "/dev/null"}" \
            --create-output-bam-index ~{create_bam} \
            --read-validation-stringency "~{validation_stringency}" \
            --duplicate-scoring-strategy "~{duplicate_scoring_strategy}" \
            --read-name-regex '~{
                if (optical_distance > 0) then read_name_regex else "null"
            }' \
            --duplicate-tagging-policy "~{tagging_policy}" \
            --optical-duplicate-pixel-distance ~{optical_distance} \
            --spark-master local[~{ncpu}]
    >>>

    output {
        File? duplicate_marked_bam = "~{prefix}.bam"
        File? duplicate_marked_bam_index = "~{prefix}.bam.bai"
        File mark_duplicates_metrics = "~{prefix}.metrics.txt"
    }

    requirements {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        maxRetries: 1
    }
}

task apply_vqsr {
    meta {
        description: "Applies variant quality score recalibration to a VCF file."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/35967650663195-ApplyVQSR"
        outputs: {
            vcf_recalibrated: "Recalibrated VCF file",
            vcf_recalibrated_index: "Index file for the recalibrated VCF",
        }
     }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index for FASTA format genome"
        reference_dict: "Dictionary file for FASTA format genome"
        vcf: "Input VCF format file on which to apply variant quality score recalibration"
        vcf_index: "VCF index file corresponding to the input VCF"
        recal_file: "Recalibration file generated by the VariantRecalibrator tool"
        tranches_file: "Tranches file generated by the VariantRecalibrator tool"
        mode: "Variant type (SNP or INDEL) to apply variant quality score recalibration for. If `BOTH`, then the same model will be applied to both SNPs and indels."
        prefix: "Prefix for the output recalibrated VCF. The extension `.recalibrated.vcf.gz` will be added."
        truth_sensitivity_filter_level: {
            description: "Truth sensitivity level at which to filter variants.",
            help: "This corresponds to the `--truth-sensitivity-filter-level` argument of GATK ApplyVQSR. Default is 99.0, which means that variants will be filtered at the threshold that achieves 99.0% sensitivity on the truth set used to build the model."
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        File vcf
        File vcf_index
        File recal_file
        File recal_file_index
        File tranches_file
        VariantMode mode = VariantMode.SNP
        String prefix = basename(vcf, ".vcf.gz")
        Float truth_sensitivity_filter_level = 99.0
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(vcf, "GB") * 2)
        + ceil(size(reference_fasta, "GB") * 2)
        + 30 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{sub(basename(reference_fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta.fa" \
            || ln -sf "~{reference_fasta}" "$ref_fasta.fa"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{reference_dict}" "$ref_fasta.dict"

        ln -sf "~{vcf}" "~{basename(vcf)}"
        ln -sf "~{vcf_index}" "~{basename(vcf_index)}"

        gatk ApplyVQSR \
            -R "$ref_fasta.fa" \
            -V "~{basename(vcf)}" \
            --recal-file "~{recal_file}" \
            --tranches-file "~{tranches_file}" \
            --truth-sensitivity-filter-level ~{truth_sensitivity_filter_level} \
            --mode "~{mode}" \
            -O "~{prefix}.recalibrated.vcf.gz"

        rm -rf "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict" "~{basename(vcf)}" "~{basename(vcf_index)}"
    >>>

    output {
        File vcf_recalibrated = "~{prefix}.recalibrated.vcf.gz"
        File vcf_recalibrated_index = "~{prefix}.recalibrated.vcf.gz.tbi"
    }

    requirements {
        cpu: 1
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        maxRetries: 1
    }
}

task variant_recalibrator {
    meta {
        description: "Generates recalibration tables for variant quality score recalibration."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/35967662766235-VariantRecalibrator"
        outputs: {
            recal_file: "Recalibration file containing the model generated by VariantRecalibrator",
            tranches_file: "Tranches file containing information about the sensitivity and specificity of the model at different score thresholds",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index for FASTA format genome"
        reference_dict: "Dictionary file for FASTA format genome"
        vcf: "Input VCF format file on which to perform variant quality score recalibration"
        vcf_index: "VCF index file corresponding to the input VCF"
        resource_vcfs: "List of VCF files corresponding to the resources in `resources`"
        resource_vcf_indices: "List of VCF index files corresponding to the VCF files in `resource_vcfs`"
        resources: "List of resources to use for building the recalibration model"
        annotations: "List of annotations to use for building the recalibration model."
        mode: "Variant type (SNP or INDEL) to build the recalibration model for. If `BOTH`, then separate models will be built for SNPs and indels."
        prefix: "Prefix for the output recalibration files. The extensions `.recal` and `.tranches` will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        File vcf
        File vcf_index
        Array[File] resource_vcfs
        Array[File] resource_vcf_indices
        Array[String] resources
        Array[String] annotations = ["QD", "MQ", "MQRankSum", "ReadPosRankSum", "FS", "SOR"]
        VariantMode mode = VariantMode.SNP
        String prefix = basename(vcf, ".vcf.gz")
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(vcf, "GB") * 2)
        + ceil(size(reference_fasta, "GB") * 2)
        + 30 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{sub(basename(reference_fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta.fa" \
            || ln -sf "~{reference_fasta}" "$ref_fasta.fa"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{reference_dict}" "$ref_fasta.dict"

        ln -sf "~{vcf}" "~{basename(vcf)}"
        ln -sf "~{vcf_index}" "~{basename(vcf_index)}"

        for vcf in ~{sep(" ", resource_vcfs)}; do
            ln -sf "~{vcf}" "$(basename "$vcf")"
        done

        for vcf_index in ~{sep(" ", resource_vcf_indices)}; do
            ln -sf "~{vcf_index}" "$(basename "$vcf_index")"
        done

        gatk VariantRecalibrator \
            -R "$ref_fasta.fa" \
            -V "~{basename(vcf)}" \
            --mode "~{mode}" \
            -O "~{prefix}.recal" \
            --tranches-file "~{prefix}.tranches" \
            --rscript-file "~{prefix}.plots.R" \
            --dont-run-rscript \
            ~{sep(" ", prefix("--resource:", resources))} \
            ~{sep(" ", prefix("-an ", quote(annotations)))}

        rm -rf "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict" "~{basename(vcf)}" "~{basename(vcf_index)}"
        for vcf in ~{sep(" ", resource_vcfs)}; do
            rm -rf "$(basename "$vcf")"
        done

        for vcf_index in ~{sep(" ", resource_vcf_indices)}; do
            rm -rf "$(basename "$vcf_index")"
        done
    >>>

    output {
        File recal_file = "~{prefix}.recal"
        File recal_index = "~{prefix}.recal.idx"
        File tranches_file = "~{prefix}.tranches"
        File rscript_file = "~{prefix}.plots.R"
    }

    requirements {
        cpu: 1
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        maxRetries: 1
    }
}

task calculate_genotype_posteriors {
    meta {
        description: "Calculates posterior genotype probabilities for a VCF file using population allele frequencies."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/35967653366299-CalculateGenotypePosteriors"
        outputs: {
            vcf_posteriors: "VCF file with updated genotype posteriors",
            vcf_posteriors_index: "Index file for the VCF with updated genotype posteriors",
        }
    }

    parameter_meta {
        vcf: "Input VCF format file on which to calculate genotype posteriors"
        vcf_index: "VCF index file corresponding to the input VCF"
        supporting_vcf: "VCF file containing population allele frequencies to use as support for calculating genotype posteriors"
        supporting_vcf_index: "VCF index file corresponding to the supporting VCF"
        prefix: "Prefix for the output VCF with updated genotype posteriors. The extension `.posteriors.vcf.gz` will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation."
    }

    input {
        File vcf
        File vcf_index
        File supporting_vcf
        File supporting_vcf_index
        String prefix = basename(vcf, ".vcf.gz")
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(vcf, "GB") * 2) + 30 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ln -sf "~{vcf}" "~{basename(vcf)}"
        ln -sf "~{vcf_index}" "~{basename(vcf_index)}"
        ln -sf "~{supporting_vcf}" "supporting.vcf.gz"
        ln -sf "~{supporting_vcf_index}" "supporting.vcf.gz.tbi"

        gatk CalculateGenotypePosteriors \
            -V "~{basename(vcf)}" \
            -supporting "supporting.vcf.gz" \
            -O "~{prefix}.posteriors.vcf.gz"

        rm -rf "~{basename(vcf)}" "~{basename(vcf_index)}" "supporting.vcf.gz" "supporting.vcf.gz.tbi"
    >>>

    output {
        File vcf_posteriors = "~{prefix}.posteriors.vcf.gz"
        File vcf_posteriors_index = "~{prefix}.posteriors.vcf.gz.tbi"
    }

    requirements {
        cpu: 1
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        maxRetries: 1
    }
}

task genotype_gvcfs {
    meta {
        description: "Generates per-sample GVCF files from BAM files using GATK's HaplotypeCaller in GVCF mode."
        external_help: "https://gatk.broadinstitute.org/hc/en-us/articles/35967678260379-GenotypeGVCFs"
        outputs: {
            vcf: "VCF file containing variant and non-variant sites with genotype likelihoods",
            vcf_index: "Index file for the VCF",
        }
    }

    parameter_meta {
        gvcf: "Input GVCF format file containing GVCF records to genotype"
        gvcf_index: "GVCF index file corresponding to the input GVCF"
        fasta: "Reference genome in FASTA format"
        fasta_index: "Index for FASTA format genome"
        dict: "Dictionary file for FASTA format genome"
        prefix: "Prefix for the output GVCF. The extension `.g.vcf.gz` will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
    }

    input {
        File gvcf
        File gvcf_index
        File fasta
        File fasta_index
        File dict
        String prefix
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(gvcf, "GB") * 2)
        + 30
        + ceil(size(fasta, "GB"))
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{sub(basename(fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{fasta}" > "$ref_fasta.fa" \
            || ln -sf "~{fasta}" "$ref_fasta.fa"
        ln -sf "~{fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{dict}" "$ref_fasta.dict"

        ln -sf "~{gvcf}" "~{basename(gvcf)}"
        ln -sf "~{gvcf_index}" "~{basename(gvcf_index)}"

        gatk \
           --java-options "-Xmx4g" \
            GenotypeGVCFs \
            -R "$ref_fasta.fa" \
            -V "~{basename(gvcf)}" \
            -O "~{prefix}.g.vcf.gz" \

        rm -rf "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict" "~{basename(gvcf)}" "~{basename(gvcf_index)}"
    >>>

    output {
        File vcf = "~{prefix}.g.vcf.gz"
        File vcf_index = "~{prefix}.g.vcf.gz.tbi"
    }

    requirements {
        cpu: 1
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
    }
}

workflow germline_variant_calling_wf {
    meta {
        description: "Workflow for calling germline variants using GATK's HaplotypeCaller."
        outputs: {
            raw_gvcf: "Raw gVCF file output by HaplotypeCaller",
            raw_gvcf_index: "Index file for the raw gVCF",
            raw_vcf: "Raw VCF file output by GenotypeGVCFs",
            raw_vcf_index: "Index file for the raw VCF",
            recalibrated_vcf: "Recalibrated VCF file after applying VQSR",
            recalibrated_vcf_index: "Index file for the recalibrated VCF",
            vcf_final: "Final VCF file after calculating genotype posteriors",
            vcf_final_index: "Index file for the final VCF",
        }
    }

    parameter_meta {
        bam: "Input BAM format file on which to perform germline variant calling"
        bam_index: "BAM index file corresponding to the input BAM"
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index for FASTA format genome"
        reference_dict: "Dictionary file for FASTA format genome"
        dbSNP_vcf: "dbSNP VCF file"
        dbSNP_vcf_index: "dbSNP VCF index file"
        interval_list: "Interval list indicating regions in which to call variants"
        known_indels_sites_vcfs: "List of VCF files containing known indel sites to use for base quality score recalibration"
        known_indels_sites_indices: "List of VCF index files corresponding to the VCF files in `known_indels_sites_vcfs`"
        resources: {
            description: "List of resources to use for building the variant quality score recalibration model.",
            help: " Each resource should be a tuple containing the name of the resource, the path to the VCF file for the resource, and the path to the VCF index file for the resource."
        }
        prefix: "Prefix for the output files."
    }

    input {
        File bam
        File bam_index
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        #@except: SnakeCase
        File dbSNP_vcf
        #@except: SnakeCase
        File dbSNP_vcf_index
        File interval_list
        Array[File] known_indels_sites_vcfs
        Array[File] known_indels_sites_indices
        Array[Resource] resources
        String prefix = basename(bam, ".bam")
    }

    scatter (resource in resources) {
        call resource_to_string {
            res = resource
        }
    }

    call base_recalibrator {
        bam,
        bam_index,
        fasta = reference_fasta,
        fasta_index = reference_fasta_index,
        dict = reference_dict,
        dbSNP_vcf,
        dbSNP_vcf_index,
        known_indels_sites_vcfs,
        known_indels_sites_indices,
        outfile_name = prefix + ".recalibration_report.txt",
    }

    call apply_bqsr {
        bam,
        bam_index,
        recalibration_report = base_recalibrator.recalibration_report,
        prefix = prefix + ".recal",
    }

    call haplotype_caller {
        bam = apply_bqsr.recalibrated_bam,
        bam_index = apply_bqsr.recalibrated_bam_index,
        interval_list,
        fasta = reference_fasta,
        fasta_index = reference_fasta_index,
        dict = reference_dict,
        dbSNP_vcf,
        dbSNP_vcf_index,
        prefix,
        reference_confidence = ref_confidence.GVCF,
    }

    call genotype_gvcfs {
        gvcf = haplotype_caller.vcf,
        gvcf_index = haplotype_caller.vcf_index,
        fasta = reference_fasta,
        fasta_index = reference_fasta_index,
        dict = reference_dict,
        prefix,
    }

     call variant_recalibrator {
        reference_fasta,
        reference_fasta_index,
        reference_dict,
        vcf = genotype_gvcfs.vcf,
        vcf_index = genotype_gvcfs.vcf_index,
        resources = resource_to_string.res_string,
        resource_vcfs = resource_to_string.res_vcf,
        resource_vcf_indices = resource_to_string.res_vcf_index,
    }

    call apply_vqsr {
        reference_fasta,
        reference_fasta_index,
        reference_dict,
        vcf = genotype_gvcfs.vcf,
        vcf_index = genotype_gvcfs.vcf_index,
        recal_file = variant_recalibrator.recal_file,
        recal_file_index = variant_recalibrator.recal_index,
        tranches_file = variant_recalibrator.tranches_file,
        prefix = prefix + ".vqsr",
    }

    call calculate_genotype_posteriors {
        vcf = apply_vqsr.vcf_recalibrated,
        vcf_index = apply_vqsr.vcf_recalibrated_index,
        supporting_vcf = dbSNP_vcf,
        supporting_vcf_index = dbSNP_vcf_index,
        prefix = prefix + ".posteriors",
    }

    output {
        File raw_gvcf = haplotype_caller.vcf
        File raw_gvcf_index = haplotype_caller.vcf_index
        File raw_vcf = genotype_gvcfs.vcf
        File raw_vcf_index = genotype_gvcfs.vcf_index
        File recalibrated_vcf = apply_vqsr.vcf_recalibrated
        File recalibrated_vcf_index = apply_vqsr.vcf_recalibrated_index
        File vcf_final = calculate_genotype_posteriors.vcf_posteriors
        File vcf_final_index = calculate_genotype_posteriors.vcf_posteriors_index
    }
}
