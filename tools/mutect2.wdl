version 1.3

workflow mutect2_wf {
    meta {
        description: "Workflow for calling somatic variants using GATK Mutect2 and filtering with FilterMutectCalls"
        outputs: {
            filtered_somatic_vcf: "VCF file with filtered somatic variants from Mutect2",
            filtered_somatic_vcf_index: "Index file for the filtered Mutect2 somatic variants VCF",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        reference_fasta_dict: "Dictionary file for the reference genome FASTA"
        normal_bam: "Input BAM file with aligned reads for normal sample"
        normal_bam_index: "Index file for the normal BAM file"
        tumor_bam: "Input BAM file with aligned reads for tumor sample"
        tumor_bam_index: "Index file for the tumor BAM file"
        variant_vcf: "VCF file with variants and allele frequencies to summarize pileups over."
        variant_vcf_index: "Index file for the variant VCF"
        intervals: "One or more genomic intervals over which to operate. Often the same file as `variants`"
        intervals_index: "Index file for the intervals file"
        germline_resource_vcf: "Optional VCF file with germline variants for Mutect2, recommended to be from gnomAD or similar population resource"
        germline_resource_vcf_index: "Index file for the germline resource VCF"
        panel_of_normals_vcf: "Optional VCF file with panel of normals for Mutect2, recommended to be generated from a large set of normal samples processed with Mutect2"
        panel_of_normals_vcf_index: "Index file for the panel of normals VCF"
        normal_sample_name: "Name of the normal sample"
        tumor_sample_name: "Name of the tumor sample"
        output_prefix: "Prefix for output files. The extensions '.vcf.gz' and '.vcf.gz.tbi' will be added."
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File reference_fasta_dict
        File normal_bam
        File normal_bam_index
        File tumor_bam
        File tumor_bam_index
        File variant_vcf
        File variant_vcf_index
        File intervals
        File intervals_index
        File? germline_resource_vcf
        File? germline_resource_vcf_index
        File? panel_of_normals_vcf
        File? panel_of_normals_vcf_index
        String normal_sample_name = basename(normal_bam, ".bam")
        String tumor_sample_name = basename(tumor_bam, ".bam")
        String output_prefix = basename(tumor_bam, ".bam") + "_v_" + basename(normal_bam, ".bam")
    }

    call mutect2 {
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
        normal_bam,
        normal_bam_index,
        tumor_bam,
        tumor_bam_index,
        normal_sample_name,
        tumor_sample_name,
        output_prefix = output_prefix + "_unfiltered",
        germline_resource_vcf,
        germline_resource_vcf_index,
        panel_of_normals_vcf,
        panel_of_normals_vcf_index,
    }

    call get_pileup_summaries as get_tumor_pileups {
        bam = tumor_bam,
        bam_index = tumor_bam_index,
        intervals,
        intervals_index,
        variants = variant_vcf,
        variants_index = variant_vcf_index,
        prefix = output_prefix + "_tumor",
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
    }

    call get_pileup_summaries as get_normal_pileups {
        bam = normal_bam,
        bam_index = normal_bam_index,
        intervals,
        intervals_index,
        variants = variant_vcf,
        variants_index = variant_vcf_index,
        prefix = output_prefix + "_normal",
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
    }

    call calculate_contamination {
        tumor_pileups = get_tumor_pileups.pileup_summaries,
        normal_pileups = get_normal_pileups.pileup_summaries,
        prefix = output_prefix,
    }

    call filter_mutect {
        unfiltered_somatic_vcf = mutect2.somatic_vcf,
        unfiltered_somatic_vcf_index = mutect2.somatic_vcf_index,
        unfiltered_somatic_vcf_stats = mutect2.stats,
        reference_fasta,
        reference_fasta_index,
        reference_fasta_dict,
        contamination_table = calculate_contamination.contamination_table,
        maf_segments = calculate_contamination.maf_segments,
        prefix = output_prefix,
    }

    output {
        File filtered_somatic_vcf = filter_mutect.filtered_somatic_vcf
        File filtered_somatic_vcf_index = filter_mutect.filtered_somatic_vcf_index
    }
}

task mutect2 {
    meta {
        description: "Run GATK Mutect2 somatic variant calling workflow"
        outputs: {
            somatic_vcf: "VCF file with somatic variants called by Mutect2",
            somatic_vcf_index: "Index file for the Mutect2 somatic variants VCF",
            stats: "VCF file with statistics about the Mutect2 somatic variant calls",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        reference_fasta_dict: "Dictionary file for the reference genome FASTA"
        normal_bam: "Input BAM file with aligned reads for normal sample"
        normal_bam_index: "Index file for the normal BAM file"
        tumor_bam: "Input BAM file with aligned reads for tumor sample"
        tumor_bam_index: "Index file for the tumor BAM file"
        germline_resource_vcf: "Optional VCF file with germline variants for Mutect2, recommended to be from gnomAD or similar population resource"
        germline_resource_vcf_index: "Index file for the germline resource VCF"
        panel_of_normals_vcf: "Optional VCF file with panel of normals for Mutect2, recommended to be generated from a large set of normal samples processed with Mutect2"
        panel_of_normals_vcf_index: "Index file for the panel of normals VCF"
        normal_sample_name: "Name of the normal sample"
        tumor_sample_name: "Name of the tumor sample"
        output_prefix: "Prefix for output file. The extension '.vcf.gz' will be added."
        threads: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File reference_fasta_dict
        File normal_bam
        File normal_bam_index
        File tumor_bam
        File tumor_bam_index
        File? germline_resource_vcf
        #@except: UnusedInput
        File? germline_resource_vcf_index
        File? panel_of_normals_vcf
        #@except: UnusedInput
        File? panel_of_normals_vcf_index
        String normal_sample_name = basename(normal_bam, ".bam")
        String tumor_sample_name = basename(tumor_bam, ".bam")
        String output_prefix = "mutect2_somatic_variants"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2)
        + ceil(size(normal_bam, "GB") * 2)
        + ceil(size(tumor_bam, "GB") * 2)
        + 20
        + modify_disk_size_gb

    String tumor = basename(tumor_bam)
    String normal = basename(normal_bam)

    command <<<
        set -euo pipefail

        ref_fasta=~{sub(basename(reference_fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta.fa" \
            || ln -sf "~{reference_fasta}" "$ref_fasta.fa"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{reference_fasta_dict}" "$ref_fasta.dict"

        ln -sf "~{tumor_bam}" "~{tumor}"
        ln -sf "~{tumor_bam_index}" "~{tumor}.bai"
        ln -sf "~{normal_bam}" "~{normal}"
        ln -sf "~{normal_bam_index}" "~{normal}.bai"

        gatk Mutect2 \
            -R "$ref_fasta.fa" \
            -I "~{tumor}" \
            -I "~{normal}" \
            -normal "~{normal_sample_name}" \
            -tumor "~{tumor_sample_name}" \
            ~{if defined(germline_resource_vcf) then "-germline-resource '" + germline_resource_vcf + "'" else ""} \
            ~{if defined(panel_of_normals_vcf) then "-panel-of-normals '" + panel_of_normals_vcf + "'" else ""} \
            -O "~{output_prefix}.vcf.gz" \
            --native-pair-hmm-threads "~{threads}"

        rm -rf "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict" "~{tumor}" "~{tumor}.bai" "~{normal}" "~{normal}.bai"
     >>>

     output {
         File somatic_vcf = "~{output_prefix}.vcf.gz"
         File somatic_vcf_index = "~{output_prefix}.vcf.gz.tbi"
         File stats = "~{output_prefix}.vcf.gz.stats"
    }

    requirements {
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        cpu: threads
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
    }
}

task filter_mutect {
    meta {
        description: "Run GATK FilterMutectCalls to filter Mutect2 somatic variant calls"
        outputs: {
            filtered_somatic_vcf: "VCF file with filtered somatic variants from Mutect2",
            filtered_somatic_vcf_index: "Index file for the filtered Mutect2 somatic variants VCF",
        }
    }

    parameter_meta {
        unfiltered_somatic_vcf: "Input VCF file with unfiltered somatic variants from Mutect2"
        unfiltered_somatic_vcf_index: "Index file for the unfiltered Mutect2 somatic variants VCF"
        unfiltered_somatic_vcf_stats: "VCF file with statistics about the unfiltered Mutect2 somatic variant calls"
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        reference_fasta_dict: "Dictionary file for the reference genome FASTA"
        contamination_table: "Table file with estimated contamination fraction and related metrics from GATK CalculateContamination"
        maf_segments: "Table file with segmented genomic intervals for MAF calculation from GATK CalculateContamination"
        prefix: "Prefix for output file. The extension '.vcf.gz' will be added."
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File unfiltered_somatic_vcf
        File unfiltered_somatic_vcf_index
        File unfiltered_somatic_vcf_stats
        File reference_fasta
        File reference_fasta_index
        File reference_fasta_dict
        File contamination_table
        File maf_segments
        String prefix = basename(unfiltered_somatic_vcf, ".vcf.gz") + "_filtered"
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2)
        + ceil(size(unfiltered_somatic_vcf, "GB") * 2)
        + 20
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{sub(basename(reference_fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta.fa" \
             || ln -sf "~{reference_fasta}" "$ref_fasta.fa"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{reference_fasta_dict}" "$ref_fasta.dict"

        ln -sf "~{unfiltered_somatic_vcf}" "~{basename(unfiltered_somatic_vcf)}"
        ln -sf "~{unfiltered_somatic_vcf_index}" "~{basename(unfiltered_somatic_vcf_index)}"
        ln -sf "~{unfiltered_somatic_vcf_stats}" "~{basename(unfiltered_somatic_vcf_stats)}"

        gatk --java-options "-Xmx~{24000}m" \
            FilterMutectCalls \
            -R "$ref_fasta.fa" \
            -V "~{basename(unfiltered_somatic_vcf)}" \
            --contamination-table "~{contamination_table}" \
            --tumor-segmentation "~{maf_segments}" \
            -O "~{prefix}.vcf.gz"

        rm -rf "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict" "~{basename(unfiltered_somatic_vcf)}" "~{basename(unfiltered_somatic_vcf_index)}" "~{basename(unfiltered_somatic_vcf_stats)}"
    >>>

    output {
        File filtered_somatic_vcf = "~{prefix}.vcf.gz"
        File filtered_somatic_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    requirements {
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        cpu: 1
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
    }
}

task calculate_contamination {
    meta {
        description: "Run GATK CalculateContamination to estimate cross-sample contamination in tumor sample"
        outputs: {
            contamination_table: "Table file with estimated contamination fraction and related metrics",
            maf_segments: "Table file with segmented genomic intervals for MAF calculation",
        }
    }

    parameter_meta {
        tumor_pileups: "Pileup summaries file for the tumor sample generated by GATK GetPileupSummaries"
        normal_pileups: {
            description: "Pileup summaries file for the normal sample generated by GATK GetPileupSummaries.",
            help: "Optional but recommended for more accurate contamination estimation."
        }
        prefix: "Prefix for output files. The extensions '.table' and '.segments.table' will be added."
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
      File tumor_pileups
      File? normal_pileups
      String prefix = basename(tumor_pileups, ".table")
      Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(tumor_pileups, "GB") * 2)
        + ceil(size(normal_pileups, "GB") * 2)
        + 20
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{24000}m" \
            CalculateContamination \
            -I "~{tumor_pileups}" \
            -O "~{prefix}.contamination.table" \
            --tumor-segmentation "~{prefix}.segments.table" \
            ~{if defined(normal_pileups) then "-matched '" + normal_pileups + "'" else ""}
    >>>

    output {
        File contamination_table = "~{prefix}.contamination.table"
        File maf_segments = "~{prefix}.segments.table"
    }

    requirements {
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        memory: "25 GB"
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
        cpu: 1
    }
}

task get_pileup_summaries {
    meta {
        description: "Run GATK GetPileupSummaries to generate pileup summaries for contamination estimation"
        outputs: {
            pileup_summaries: "Table file with pileup summaries for the input BAM file over the specified variants and intervals"
        }
    }

    parameter_meta {
        bam: "Input BAM file with aligned reads"
        bam_index: "Index file for the input BAM file"
        intervals: "One or more genomic intervals over which to operate. Often the same file as `variants`"
        intervals_index: "Index file for the intervals file"
        variants: "VCF file with variants and allele frequencies to summarize pileups over."
        variants_index: "Index file for the variant VCF"
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        reference_fasta_dict: "Dictionary file for the reference genome FASTA"
        prefix: "Prefix for output file. The extension '.table' will be added."
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File bam
        File bam_index
        File intervals
        File intervals_index
        File variants
        File variants_index
        File reference_fasta
        File reference_fasta_index
        File reference_fasta_dict
        String prefix = basename(bam, ".bam") + "_pileup_summaries"
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(bam, "GB") * 2)
        + ceil(size(intervals, "GB") * 2)
        + ceil(size(variants, "GB") * 2)
        + ceil(size(reference_fasta, "GB") * 2)
        + 20
        + modify_disk_size_gb

    String name = basename(bam, ".bam")

    command <<<
        set -euo pipefail

        ln -sf "~{bam}" "~{name}.bam"
        ln -sf "~{bam_index}" "~{name}.bam.bai"
        ln -sf "~{intervals}" "~{basename(intervals)}"
        ln -sf "~{intervals_index}" "~{basename(intervals_index)}"
        ln -sf "~{variants}" "~{basename(variants)}"
        ln -sf "~{variants_index}" "~{basename(variants_index)}"
        ref_fasta=~{sub(basename(reference_fasta, ".gz"), ".(fasta|fa)?$", "")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta.fa" \
            || ln -sf "~{reference_fasta}" "$ref_fasta.fa"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fa.fai"
        ln -sf "~{reference_fasta_dict}" "$ref_fasta.dict"

        gatk --java-options "-Xmx~{(task.memory / (1024 * 1024)) - 3000}m" \
            GetPileupSummaries \
            --disable-bam-index-caching \
            -I "~{name}.bam" \
            -V "~{basename(variants)}" \
            -L "~{basename(intervals)}" \
            -O "~{prefix}.table"

        rm -rf "~{name}.bam" "~{name}.bam.bai" "~{basename(intervals)}" "~{basename(intervals_index)}" "~{basename(variants)}" "~{basename(variants_index)}" "$ref_fasta.fa" "$ref_fasta.fa.fai" "$ref_fasta.dict"
    >>>

    output {
        File pileup_summaries = "~{prefix}.table"
    }

    requirements {
        container: "quay.io/biocontainers/gatk4:4.6.2.0--py310hdfd78af_1"
        memory: "50 GB"
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
        cpu: 1
    }
}
