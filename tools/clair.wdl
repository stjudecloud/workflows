version 1.3

# [ClairS pre-trained models](https://github.com/HKU-BAL/ClairS#pre-trained-models)
#
# | `--platform`                   | Model                     | Chemistry / Instrument              | Basecaller                                                |
# | ------------------------------ | ------------------------- | ----------------------------------- | --------------------------------------------------------- |
# | `ont_r10_dorado_sup_4khz`      | r1041_e82_400bps_sup_v410 | R10.4.1, 4kHz                       | Dorado SUP                                                |
# | `ont_r10_dorado_sup_5khz_ss`   | r1041_e82_400bps_sup_v420 | R10.4.1, 5kHz                       | Dorado SUP (synthetic-sample trained)                     |
# | `ont_r10_dorado_sup_5khz_ssrs` | r1041_e82_400bps_sup_v420 | R10.4.1, 5kHz                       | Dorado SUP (synthetic + real-sample trained; recommended) |
# | `ont_r10_dorado_hac_5khz`      | r1041_e82_400bps_hac_v420 | R10.4.1, 5kHz                       | Dorado HAC                                                |
# | `ont_r10_dorado_hac_4khz`      | r1041_e82_400bps_hac_v410 | R10.4.1, 4kHz                       | Dorado HAC                                                |
# | `ont_r10_guppy`                | r104_e81_sup_g5015        | R10.4/R10.4.1, 4kHz                 | Guppy5 SUP                                                |
# | `ont_r9_guppy`                 | r941_prom_sup_g5014       | R9.4.1, 4kHz                        | Guppy5 SUP                                                |
# | `ilmn`                         | ilmn                      | NovaSeq / HiSeq X                   | n/a                                                       |
# | `hifi_sequel2`                 | hifi_sequel2              | PacBio Sequel II, Chemistry 2.0     | n/a                                                       |
# | `hifi_revio`                   | hifi_revio                | PacBio Revio, SMRTbell prep kit 3.0 | n/a                                                       |
enum ClairSPlatform {
    ont_r10_dorado_sup_4khz,
    ont_r10_dorado_sup_5khz_ss,
    ont_r10_dorado_sup_5khz_ssrs,
    ont_r10_dorado_hac_5khz,
    ont_r10_dorado_hac_4khz,
    ont_r10_guppy,
    ont_r9_guppy,
    ilmn,
    hifi_sequel2,
    hifi_revio,
}

enum Clair3Platform {
    ont,
    hifi,
    ilmn,
}

# [Clair3 pre-trained models](https://github.com/HKU-BAL/Clair3#pre-trained-models)
#
# | Model                                          | Platform | Chemistry / Basecaller                 | Notes                                        |
# | ---------------------------------------------- | -------- | -------------------------------------- | -------------------------------------------- |
# | `ont`                                          | ont      | ONT (legacy/generic default)           |                                              |
# | `ont_guppy5`                                   | ont      | R9.4.1, Guppy5 SUP                     | Also usable on HAC reads                     |
# | `r941_prom_sup_g5014`                          | ont      | R9.4.1, Guppy5 SUP                     |                                              |
# | `r941_prom_hac_g360+g422`                      | ont      | R9.4.1, Guppy3/4 HAC                   |                                              |
# | `r1041_e82_400bps_hac_v410`                    | ont      | R10.4.1 E8.2 (4kHz), Dorado v4.1.0 HAC |                                              |
# | `r1041_e82_400bps_hac_v500`                    | ont      | R10.4.1 E8.2 (5kHz), Dorado v5.0.0 HAC |                                              |
# | `r1041_e82_400bps_hac_v520`                    | ont      | R10.4.1 E8.2 (5kHz), Dorado v5.2.0 HAC |                                              |
# | `r1041_e82_400bps_hac_v520_with_mv`            | ont      | R10.4.1 E8.2 (5kHz), Dorado v5.2.0 HAC | Signal-aware; requires `--enable_dwell_time` |
# | `r1041_e82_400bps_hac_v600`                    | ont      | R10.4.1 E8.2 (5kHz), Dorado v6.0.0 HAC | Latest HAC model                             |
# | `r1041_e82_400bps_hac_v600_with_mv`            | ont      | R10.4.1 E8.2 (5kHz), Dorado v6.0.0 HAC | Signal-aware; requires `--enable_dwell_time` |
# | `r1041_e82_400bps_hac_with_mv`                 | ont      | R10.4.1 E8.2, HAC                      | Signal-aware; requires `--enable_dwell_time` |
# | `r1041_e82_400bps_sup_v410`                    | ont      | R10.4.1 E8.2 (4kHz), Dorado v4.1.0 SUP |                                              |
# | `r1041_e82_400bps_sup_v430_bacteria_finetuned` | ont      | R10.4.1 E8.2, Dorado v4.3.0 SUP        | Fine-tuned on 12 bacterial genomes           |
# | `r1041_e82_400bps_sup_v500`                    | ont      | R10.4.1 E8.2 (5kHz), Dorado v5.0.0 SUP |                                              |
# | `r1041_e82_400bps_sup_v520`                    | ont      | R10.4.1 E8.2 (5kHz), Dorado v5.2.0 SUP | Latest SUP model                             |
# | `r1041_e82_400bps_sup_v520_with_mv`            | ont      | R10.4.1 E8.2 (5kHz), Dorado v5.2.0 SUP | Signal-aware; requires `--enable_dwell_time` |
# | `r1041_e82_400bps_sup_with_mv`                 | ont      | R10.4.1 E8.2, SUP                      | Signal-aware; requires `--enable_dwell_time` |
# | `hifi`                                         | hifi     | PacBio HiFi (legacy/generic default)   |                                              |
# | `hifi_sequel2`                                 | hifi     | PacBio HiFi, Sequel II                 |                                              |
# | `hifi_revio`                                   | hifi     | PacBio HiFi, Revio                     |                                              |
# | `ilmn`                                         | ilmn     | Illumina short-read                    |                                              |
enum Clair3Model[String] {
    hifi,
    hifi_sequel2,
    ont,
    r1041_e82_400bps_hac_v410,
    r1041_e82_400bps_hac_v500,
    r1041_e82_400bps_hac_v520,
    r1041_e82_400bps_hac_v520_with_mv,
    r1041_e82_400bps_hac_v600,
    r1041_e82_400bps_hac_v600_with_mv,
    r1041_e82_400bps_hac_with_mv,
    r1041_e82_400bps_sup_v410,
    r1041_e82_400bps_sup_v430_bacteria_finetuned,
    r1041_e82_400bps_sup_v500,
    r1041_e82_400bps_sup_v520,
    r1041_e82_400bps_sup_v520_with_mv,
    r1041_e82_400bps_sup_with_mv,
    r941_prom_sup_g5014,
    hifi_revio,
    ilmn,
    ont_guppy5,
    r941_prom_hac_g360_g422 = "r941_prom_hac_g360+g422",
}

task clair3 {
    meta {
        description: "Run Clair3 variant caller for small variants using deep neural networks"
        outputs: {
            pileup_vcf: "VCF file with variants called using pileup model",
            full_alignment_vcf: "VCF file with variants called using full-alignment model. Undefined if no full-alignment candidates were selected (e.g. too few or no low-quality variants in the pileup output for the called region).",
            merged_vcf: "Final merged VCF file with variants from both models",
        }
    }

    parameter_meta {
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        bam: "Input BAM file with aligned reads"
        bam_index: "Index file for the input BAM file"
        model: {
            description: "Pre-trained Clair3 model to use for variant calling",
        }
        bed_regions: "Optional BED file specifying regions to call variants in"
        vcf_candidates: "Optional VCF file with candidate variants to consider"
        contigs: "Optional list of contigs to call variants in. If undefined, all contigs in the reference FASTA will be considered."
        output_dir: "Directory to store Clair3 output"
        platform: "Sequencing platform used to generate the reads"
        all_contigs: "Boolean indicating whether to include all contigs in variant calling. If false only chr{1..22,X,Y} are called."
        print_ref_calls: "Boolean indicating whether to print reference calls in the output VCF"
        gvcf: "Boolean indicating whether to output gVCF format"
        ncpu: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
    }

    input {
        File reference_fasta
        File reference_fasta_index
        File bam
        File bam_index
        File? bed_regions
        File? vcf_candidates
        Array[String] contigs = []
        Clair3Model model = Clair3Model.ilmn
        Clair3Platform platform = Clair3Platform.ilmn
        String output_dir = "clair3_output"
        Boolean all_contigs = false
        Boolean print_ref_calls = false
        Boolean gvcf = false
        Int ncpu = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2) + ceil(size(bam, "GB") * 2) + 150
        + modify_disk_size_gb

    String filename = basename(bam)

    #@ except: ShellCheck
    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        # Clair-3 resolves the BAM path...
        ln -s  "~{bam}" "~{filename}"
        ln -s  "~{bam_index}" "~{filename}.bai"

        run_clair3.sh \
            --bam_fn="~{filename}" \
            --ref_fn="$ref_fasta" \
            --threads="~{ncpu}" \
            --platform="~{platform}" \
            --model_path="/opt/models/~{model}" \
            --output="~{output_dir}" \
            ~{if length(contigs) > 0 then "--ctg_name='~{sep(",", contigs)}'" else ""} \
            ~{if all_contigs then "--include_all_ctgs" else ""} \
            ~{if print_ref_calls then "--print_ref_calls" else ""} \
            ~{if defined(bed_regions) then "--bed_fn='~{bed_regions}'" else ""} \
            ~{if defined(vcf_candidates) then "--vcf_fn='~{vcf_candidates}'" else ""} \
            ~{if gvcf then "--gvcf" else ""}

        rm -rf "$ref_fasta" "$ref_fasta.fai" "~{filename}" "~{filename}.bai"
    >>>

    output {
        File pileup_vcf = "~{output_dir}/pileup.vcf.gz"
        File? full_alignment_vcf = "~{output_dir}/full_alignment.vcf.gz"
        File merged_vcf = "~{output_dir}/merge_output.vcf.gz"
    }

    requirements {
        container: "hkubal/clair3:v2.0.3"
        cpu: ncpu
        memory: "64 GB"
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}

task clair_s {
    meta {
        description: "Run ClairS paired sample variant caller"
        outputs: {
            vcf: "VCF file with somatic variants called by ClairS",
            vcf_index: "Index for `vcf`",
        }
    }

    parameter_meta {
        tumor_bam: "Input BAM file with aligned reads for tumor sample"
        tumor_bam_index: "Index file for the tumor BAM file"
        normal_bam: "Input BAM file with aligned reads for normal sample"
        normal_bam_index: "Index file for the normal BAM file"
        reference_fasta: "Reference genome in FASTA format"
        reference_fasta_index: "Index file for the reference genome FASTA"
        platform: "Sequencing platform used to generate the reads"
        bed_regions: "Optional BED file specifying regions to call variants in"
        vcf_candidates: "Optional VCF file with candidate variants to consider"
        pileup_model: "Optional pre-trained ClairS pileup model to use for variant calling"
        full_alignment_model: "Optional pre-trained ClairS full-alignment model to use for variant calling"
        snv_min_qual: "Minimum quality score required to call a SNV"
        indel_min_qual: "Minimum quality score required to call an indel"
        contigs: "Optional list of contigs to call variants in. If undefined, all contigs in the reference FASTA will be considered."
        prefix: "Prefix for ClairS output files"
        sample_name: "Sample name to use in the output VCF"
        output_dir: "Directory to store ClairS output"
        all_contigs: "Boolean indicating whether to include all contigs in variant calling. If false only chr{1..22,X,Y} are called."
        print_ref_calls: "Boolean indicating whether to print reference calls in the output VCF"
        print_germline_calls: "Boolean indicating whether to print germline calls in the output VCF"
        remove_intermediate_directory: "Boolean indicating whether to remove the intermediate output directory generated by ClairS after the final VCF is produced"
        snv_min_af: "Minimum allele frequency required to call a SNV"
        ncpu: "Number of threads to use"
        modify_disk_size_gb: "Additional disk size in GB to allocate"
        min_coverage: "Minimum coverage required at a site to make a variant call"
        chunk_size: "Size of genomic chunks in base pairs to split the input BAMs into for parallel processing"
    }

    input {
        File tumor_bam
        File tumor_bam_index
        File normal_bam
        File normal_bam_index
        File reference_fasta
        File reference_fasta_index
        ClairSPlatform platform
        File? bed_regions
        File? vcf_candidates
        File? pileup_model
        File? full_alignment_model
        Int? snv_min_qual
        Int? indel_min_qual
        Array[String] contigs = []
        String prefix = "output"
        String sample_name = "SAMPLE"
        String output_dir = "output"
        Boolean all_contigs = false
        Boolean print_ref_calls = false
        Boolean print_germline_calls = false
        Boolean remove_intermediate_directory = false
        Float snv_min_af = 0.05
        Int ncpu = 4
        Int modify_disk_size_gb = 0
        Int min_coverage = 4
        Int chunk_size = 5000000
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GB") * 2) + ceil(size(tumor_bam, "GB")
        * 2) + ceil(size(normal_bam, "GB") * 2) + 50 + modify_disk_size_gb

    String tumor = basename(tumor_bam)
    String normal = basename(normal_bam)

    #@ except: ShellCheck
    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"
        ln -sf "~{reference_fasta_index}" "$ref_fasta.fai"

        ln -s "~{tumor_bam}" "~{tumor}"
        ln -s "~{tumor_bam_index}" "~{tumor}.bai"
        ln -s "~{normal_bam}" "~{normal}"
        ln -s "~{normal_bam_index}" "~{normal}.bai"

        #@except: ShellCheck
        run_clairs \
            --tumor_bam_fn "~{tumor}" \
            --normal_bam_fn "~{normal}" \
            --ref_fn "$ref_fasta" \
            --threads "~{ncpu}" \
            --platform "~{platform}" \
            --output_dir "~{output_dir}" \
            --output_prefix "~{prefix}" \
            --min_coverage "~{min_coverage}" \
            --chunk_size "~{chunk_size}" \
            --snv_min_af "~{snv_min_af}" \
            --sample_name "~{sample_name}" \
            ~{if length(contigs) > 0 then "--ctg_name='~{sep(",", contigs)}'" else ""} \
            ~{if defined(bed_regions) then "--bed_fn '~{bed_regions}'" else ""} \
            ~{if defined(vcf_candidates)
                then "--genotyping_mode_vcf_fn '~{vcf_candidates}'"
                else ""} \
            ~{if all_contigs then "--include_all_ctgs" else ""} \
            ~{if print_ref_calls then "--print_ref_calls" else ""} \
            ~{if print_germline_calls then "--print_germline_calls" else ""} \
            ~{if defined(pileup_model)
                then "--pileup_model_path '~{pileup_model}'"
                else ""} \
            ~{if defined(full_alignment_model)
                then "--full_alignment_model_path '~{full_alignment_model}'"
                else ""} \
            ~{if defined(snv_min_qual) then "--snv_min_qual '~{snv_min_qual}'" else ""} \
            ~{if defined(indel_min_qual)
                then "--indel_min_qual '~{indel_min_qual}'"
                else ""} \
            ~{if remove_intermediate_directory
                then "--remove_intermediate_directory"
                else ""}

        rm -rf "$ref_fasta" "$ref_fasta.fai" \
            "~{tumor}" "~{tumor}.bai" \
            "~{normal}" "~{normal}.bai"
    >>>

    output {
        File vcf = "~{output_dir}/~{prefix}.vcf.gz"
        File vcf_index = "~{output_dir}/~{prefix}.vcf.gz.tbi"
    }

    requirements {
        container: "hkubal/clairs:v0.5.1"
        cpu: ncpu
        memory: "64 GB"
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}
