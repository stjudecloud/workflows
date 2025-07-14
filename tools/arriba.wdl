## [Homepage](https://arriba.readthedocs.io/en/latest/)

version 1.1

task arriba {
    meta {
        description: "Run Arriba structural variant caller on a RNA-Seq BAM file."
        help: "Typical input is a STAR-aligned BAM. Arriba also supports DRAGEN-aligned BAMs and any spec compliant BAM. That is discordant mates must have `BAM_FPROPER_PAIR (0x2)`, split reads must have `BAM_FSUPPLEMENTARY (0x800)`, and the anchor read must have a `SA` tag. Arriba also uses the `HI` tag to group supplementary alignments."
        outputs: {
            fusions: "Output file of fusions in TSV format",
            discarded_fusions: "Output file of discarded fusions in TSV format",
        }
    }

    parameter_meta {
        bam: "Input BAM format file from which to call fusions"
        gtf: "GTF features file. Gzipped or uncompressed."
        reference_fasta_gz: "Gzipped reference genome in FASTA format"
        chimeric_sam: "Optional input file of chimeric reads in SAM format, from older versions of STAR"
        exclude_list: {
            description: "Optional input file of regions to exclude from analysis in tab delimited format",
            external_help: "https://arriba.readthedocs.io/en/v2.4.0/input-files/#blacklist",
        }
        known_fusions: {
            description: "Optional input file of known fusions in tab delimited format",
            external_help: "https://arriba.readthedocs.io/en/v2.4.0/input-files/#known-fusions",
        }
        annotate_fusions: {
            description: "Optional input file in tab delimited format of fusions to annotate with tags",
            external_help: "https://arriba.readthedocs.io/en/v2.4.0/input-files/#tags",
        }
        protein_domains: {
            description: "Optional input file of protein domains coordinates in GFF3 format",
            external_help: "https://arriba.readthedocs.io/en/v2.4.0/input-files/#protein-domains",
        }
        wgs_svs: {
            description: "Optional input file of structural variants found by WGS in tab delimited or VCF format",
            external_help: "https://arriba.readthedocs.io/en/v2.4.0/input-files/#structural-variant-calls-from-wgs",
        }
        interesting_contigs: "Array of contigs to consider for analysis. Contigs can be specified with or without the prefix `chr`."
        viral_contigs: "Array of contigs to consider for viral integration site analysis."
        disable_filters: {
            description: "Array of filters to disable.",
            choices: [
                "top_expressed_viral_contigs",
                "viral_contigs",
                "low_coverage_viral_contigs",
                "uninteresting_contigs",
                "no_genomic_support",
                "short_anchor",
                "select_best",
                "many_spliced",
                "long_gap",
                "merge_adjacent",
                "hairpin",
                "small_insert_size",
                "same_gene",
                "genomic_support",
                "read_through",
                "no_coverage",
                "mismatches",
                "homopolymer",
                "low_entropy",
                "multimappers",
                "inconsistently_clipped",
                "duplicates",
                "homologs",
                "blacklist",
                "mismappers",
                "spliced",
                "relative_support",
                "min_support",
                "known_fusions",
                "end_to_end",
                "non_coding_neighbors",
                "isoforms",
                "intronic",
                "in_vitro",
                "intragenic_exonic",
                "internal_tandem_duplication",
            ],
        }
        feature_name: {
            description: "List of feature names to use in GTF.",
            help: "The Arriba default it designed to handle RefSeq, GENCODE, or ENSEMBL format annotations. `feature_name` expects a string of space/comma separated options. The required fields are `gene_name`, `gene_id`, `transcript_id`, `feature_exon`, and `feature_CDS`. The fields should space separated. The values should be provided with `field=value`. Mutliple values can be provided and separated by a pipe (`|`), e.g. `=value1|value2`. A complete example is `gene_name=gene_name|gene_id gene_id=gene_id transcript_id=transcript_id feature_exon=exon feature_CDS=CDS`.",
            external_help: "https://arriba.readthedocs.io/en/v2.4.0/command-line-options/",
            common: false,
        }
        prefix: "Prefix for the fusion result files. The extensions `.tsv` and `.discarded.tsv` will be added."
        strandedness: {
            description: "Strandedness of the input data.",
            external_help: "https://arriba.readthedocs.io/en/v2.4.0/command-line-options/",
            choices: [
                "auto",
                "yes",
                "no",
                "reverse",
            ],
        }
        mark_duplicates: {
            description: "Mark duplicates in the input BAM file with Arriba.",
            help: "Arriba performs marking of duplicates internally based on identical mapping coordinates. When this switch is set, internal marking of duplicates is disabled and Arriba assumes that duplicates have been marked by a preceding program. In this case, Arriba only discards alignments flagged with the BAM_FDUP flag. This makes sense when duplicates cannot be reliably identified solely based on their mapping coordinates, e.g. when unique molecular identifiers (UMIs) are used or when independently generated libraries are merged in a single BAM file and the read group must be interrogated to distinguish duplicates from reads that map to the same coordinates by chance. In addition, when this switch is set, duplicate reads are not considered for the calculation of the coverage at fusion breakpoints (columns coverage1 and coverage2 in the output file).",
        }
        report_additional_columns: "Report additional columns ['fusion_transcript', 'peptide_sequence', 'read_identifiers'] in the discarded fusions file."
        fill_gaps: "Fill gaps in assembled transcripts with reference bases. Expands the fusion sequence to the complete sequence of the fusion gene."
        max_e_value: "Maximum E-value for read support."
        max_mismappers: "Maximum fraction of mismapped reads in the fusion region."
        max_homolog_identity: "Maximum fraction of homologous sequence for genes."
        max_kmer_content: "Maximum fraction of repetitive 3-mer content in the fusion region."
        max_mismatch_pvalue: "Maximum p-value for mismatches in the fusion region."
        quantile: "Genes with expression above the given quantile are eligible for filtering."
        exonic_fraction: "Minimum fraction of exonic sequence between breakpoints."
        coverage_fraction: "Minimum fraction of viral contig transcription."
        min_itd_allele_fraction: "Minimum supporting read fraction for internal tandem duplications."
        max_genomic_breakpoint_distance: "With 'wgs_svs', threshold for relating genomic and transcriptomic events."
        min_supporting_reads: "Minimum number of supporting reads for a fusion."
        homopolymer_length: "Maximum homopolymer length adjacent to breakpoints."
        read_through_distance: "Minimum distance between breakpoints for read-through events."
        min_anchor_length: "Minimum anchor length for split reads."
        many_spliced_events: "Recover fusions with at least this many spliced breakpoints."
        fragment_length: "For single-end data, this is the fragment length. With paired-end reads, this is ignored and determined automatically."
        max_reads: "Subsample fusions with more than this number of reads."
        top_n: "Only report the top N most highly expressed viral integration sites."
        max_itd_length: "Maximum length of internal tandem duplications."
        min_itd_supporting_reads: "Minimum number of supporting reads for internal tandem duplications."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File gtf
        File reference_fasta_gz
        File? chimeric_sam
        File? exclude_list
        File? known_fusions
        File? annotate_fusions
        File? protein_domains
        File? wgs_svs
        Array[String] interesting_contigs = [
            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
            "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "AC_*", "NC_*",
        ]
        Array[String] viral_contigs = ["AC_*", "NC_*"]
        Array[String] disable_filters = []
        #@ except: LineWidth
        String feature_name
            = "gene_name=gene_name|gene_id,gene_id=gene_id,transcript_id=transcript_id,feature_exon=exon,feature_CDS=CDS"
        String prefix = basename(bam, ".bam") + ".fusions"
        String strandedness = "auto"
        Boolean mark_duplicates = true
        Boolean report_additional_columns = false
        Boolean fill_gaps = false
        Float max_e_value = 0.3
        Float max_mismappers = 0.8
        Float max_homolog_identity = 0.3
        Float max_kmer_content = 0.6
        Float max_mismatch_pvalue = 0.01
        Float quantile = 0.998
        Float exonic_fraction = 0.33
        Float coverage_fraction = 0.05
        Float min_itd_allele_fraction = 0.07
        Int max_genomic_breakpoint_distance = 1000000
        Int min_supporting_reads = 2
        Int homopolymer_length = 6
        Int read_through_distance = 10000
        Int min_anchor_length = 23
        Int many_spliced_events = 4
        Int fragment_length = 200
        Int max_reads = 300
        Int top_n = 5
        Int max_itd_length = 100
        Int min_itd_supporting_reads = 10
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Int bam_size_gb = ceil(size(bam, "GiB"))
    Int disk_size_gb = bam_size_gb
        + ceil(size(gtf, "GiB"))
        + ceil(size(reference_fasta_gz, "GiB"))
        + modify_disk_size_gb
    Int memory_gb = bam_size_gb + modify_memory_gb

    command <<<
        arriba \
            -x "~{bam}" \
            ~{"-c '" + chimeric_sam + "'"} \
            -o "~{prefix}.tsv" \
            -O "~{prefix}.discarded.tsv" \
            -a "~{reference_fasta_gz}" \
            -g "~{gtf}" \
            -G "~{feature_name}" \
            ~{"-b '" + exclude_list + "'"} \
            ~{"-k '" + known_fusions + "'"} \
            ~{"-t '" + annotate_fusions + "'"} \
            ~{"-p '" + protein_domains + "'"} \
            ~{"-d '" + wgs_svs + "'"} \
            -D ~{max_genomic_breakpoint_distance} \
            -s "~{strandedness}" \
            ~{(
                if length(interesting_contigs) > 0
                then "-i " + sep(",", quote(interesting_contigs))
                else ""
            )} \
            ~{(
                if length(viral_contigs) > 0
                then "-v " + sep(",", quote(viral_contigs))
                else ""
            )} \
            ~{(
                if length(disable_filters) > 0
                then "-f " + sep(",", quote(disable_filters))
                else ""
            )} \
            -E ~{max_e_value} \
            -S ~{min_supporting_reads} \
            -m ~{max_mismappers} \
            -L ~{max_homolog_identity} \
            -H ~{homopolymer_length} \
            -R ~{read_through_distance} \
            -A ~{min_anchor_length} \
            -M ~{many_spliced_events} \
            -K ~{max_kmer_content} \
            -V ~{max_mismatch_pvalue} \
            -F ~{fragment_length} \
            -U ~{max_reads} \
            -Q ~{quantile} \
            -e ~{exonic_fraction} \
            -T ~{top_n} \
            -C ~{coverage_fraction} \
            -l ~{max_itd_length} \
            -z ~{min_itd_allele_fraction} \
            -Z ~{min_itd_supporting_reads} \
            ~{if mark_duplicates then "" else "-u"} \
            ~{if report_additional_columns then "-X" else ""} \
            ~{if fill_gaps then "-I" else ""}
    >>>

    output {
        File fusions = "~{prefix}.tsv"
        File discarded_fusions = "~{prefix}.discarded.tsv"
    }

    runtime {
        cpu: 1
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/arriba:2.4.0--h0033a41_2"
        maxRetries: 1
    }
}

task arriba_tsv_to_vcf {
    meta {
        description: "Convert Arriba TSV format fusions to VCF format."
        outputs: {
            fusions_vcf: "Output file of fusions in VCF format"
        }
    }

    parameter_meta {
        fusions: "Input fusions in TSV format to convert to VCF"
        reference_fasta: "Reference genome in FASTA format. Either gzipped or uncompressed."
        prefix: "Output file name for fusions in VCF format. The extension `.vcf` will be appended."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File fusions
        File reference_fasta
        String prefix = basename(fusions, ".tsv")
        Int modify_disk_size_gb = 0
    }

    Int input_size_gb = ceil(size(fusions, "GiB"))
    Int disk_size_gb = ceil(input_size_gb)
        + (ceil(size(reference_fasta, "GiB")) * 3)
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        fasta_name=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$fasta_name" \
            || ln -sf "~{reference_fasta}" "$fasta_name"

        convert_fusions_to_vcf.sh \
            "$fasta_name" \
            "~{fusions}" \
            "~{prefix}.vcf"
    >>>

    output {
        File fusions_vcf = "~{prefix}.vcf"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/arriba:2.4.0--h0033a41_2"
        maxRetries: 1
    }
}

task arriba_extract_fusion_supporting_alignments {
    meta {
        description: "Extract alignments that support fusions."
        outputs: {
            fusion_bams: "Array of BAM files corresponding with fusions in the input file",
            fusion_bam_indexes: "Array of BAM indexes corresponding with the BAMs in the 'fusion_bams'",
        }
    }

    parameter_meta {
        bam: "Input BAM format file from which fusions were called"
        bam_index: "BAM index file corresponding to the input BAM"
        fusions: "Input fusions in TSV format for which to extract supporting alignments"
        prefix: "Output file name prefix for the extracted BAM files. The extension `.bam` will be appended."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        File bam_index
        File fusions
        String prefix = basename(fusions, ".tsv")
        Int modify_disk_size_gb = 0
    }

    Int input_size_gb = ceil(size(bam, "GiB"))
    Int disk_size_gb = ceil(input_size_gb) + 5 + modify_disk_size_gb

    command <<<
        extract_fusion-supporting_alignments.sh \
            "~{fusions}" \
            "~{bam}" \
            "~{prefix}"
    >>>

    output {
        Array[File] fusion_bams = glob("~{prefix}_*.bam")
        Array[File] fusion_bam_indexes = glob("~{prefix}_*.bam.bai")
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/arriba:2.4.0--h0033a41_2"
        maxRetries: 1
    }
}

task arriba_annotate_exon_numbers {
    meta {
        description: "Annotate fusions with exon numbers."
        outputs: {
            fusion_tsv: "TSV file with fusions annotated with exon numbers"
        }
    }

    parameter_meta {
        fusions: "Input fusions in TSV format for which to annotate gene exon numbers"
        gtf: "GTF features file. Gzipped or uncompressed."
        prefix: "Output file name for annotated fusions in TSV format. The extension `.annotated.tsv` will be appended."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File fusions
        File gtf
        String prefix = basename(fusions, ".tsv")
        Int modify_disk_size_gb = 0
    }

    Int input_size_gb = ceil(size(gtf, "GiB"))
    Int disk_size_gb = ceil(input_size_gb) + 5 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        gtf_name=~{basename(gtf, ".gz")}
        gunzip -c "~{gtf}" > "$gtf_name" || ln -sf "~{gtf}" "$gtf_name"

        annotate_exon_numbers.sh \
            "~{fusions}" \
            "$gtf_name" \
            "~{prefix}.annotated.tsv"
    >>>

    output {
        File fusion_tsv = "~{prefix}.annotated.tsv"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/arriba:2.4.0--h0033a41_2"
        maxRetries: 1
    }
}
