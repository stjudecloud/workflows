version 1.2

task align {
    meta {
        description: "Align RNA-seq reads against a reference genome using HISAT2"
        outputs: {
            alignments: "The output alignment file in SAM format"
        }
    }

    parameter_meta {
        read_one_fastq_gz: "Input gzipped FASTQ read one file to align with HISAT2"
        reference_index: "The HISAT2 index files for the reference genome"
        read_two_fastq_gz: "Input gzipped FASTQ read two file to align with HISAT2"
        output_name: "The name of the output alignment file"
        threads: "Number of threads to use for alignment"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File read_one_fastq_gz
        File reference_index
        File? read_two_fastq_gz
        String output_name = "aligned.sam"
        Int threads = 4
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil((
            size(read_one_fastq_gz, "GiB") + size(read_two_fastq_gz, "GiB")
        ) * 2)
        + ceil(size(reference_index, "GiB") * 5)
        + 10
        + modify_disk_size_gb

    command <<<
        set -euo pipefail

        mkdir hisat2_db
        tar -C hisat2_db -xzf "~{reference_index}" --no-same-owner
        PREFIX=$(basename hisat2_db/*.1.ht2 ".1.ht2")

        hisat2 \
            -q \
            -p ~{threads} \
            -S "~{output_name}" \
            -x "hisat2_db/$PREFIX" \
            -1 "~{read_one_fastq_gz}" \
            ~{if defined(read_two_fastq_gz) then "-2 \"~{read_two_fastq_gz}\"" else ""}

        rm -r hisat2_db
    >>>

    output {
        File alignments = "~{output_name}"
    }

    requirements {
        cpu: threads
        memory: "64 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/hisat2:branch-minimap2-2.2.1-0"
    }
}

task index {
    meta {
        description: "Index a reference genome for alignment with HISAT2"
        outputs: {
            reference_index: "The HISAT2 index files for the reference genome"
        }
    }

    parameter_meta {
        reference_fasta: "The reference genome in FASTA format to be indexed"
        snp: "List of SNPs"
        haplotype: "List of haplotypes"
        splice_site: "List of splice sites. Use with `exon`."
        exon: "List of exons. Use with `splice_site`."
        repeat_ref: ""
        repeat_info: ""
        repeat_snp: ""
        repeat_haplotype: ""
        bmax: "Maximum number of suffixes allowed in a block"
        seed: "Seed for psuedo-random number generator"
        bmaxdivn: "Maximum number of suffixes allowed in a block, expressed as a fraction of the length of the reference"
        index_base_name: "The base name for the output index files"
        force_large_index: "Force creation of a large index"
        disable_auto_fitting: "Disable automatic fitting of index parameters"
        nodc: "Disable difference-cover sample"
        no_ref: "Do not build bitpacked version of reference sequence for paired-end alignment"
        just_ref: "Build only the bitpacked version of reference sequence for paired-end alignment"
        threads: "Number of threads to use for indexing"
        dcv: "Period for the difference-cover sample. A larger period uses less memory, but may be slower. Must be a power of 2, no greater than 4096."
        offrate: "The off-rate for the FM index"
        ftabchars: "The lookup table to calculate initial BW range with respect to the first N characters of the query"
        localoffrate: "The off-rate for the local FM index"
        localftabchars: "The lookup table to calculate initial BW range for the local FM"
        modify_disk_size_gb: "Additional disk space to allocate (in GB)"
    }

    input {
        File reference_fasta
        File? snp
        File? haplotype
        File? splice_site
        File? exon
        File? repeat_ref
        File? repeat_info
        File? repeat_snp
        File? repeat_haplotype
        Int? bmax
        Int? seed
        Int? bmaxdivn = 4
        String index_base_name = "hisat2_index"
        Boolean force_large_index = false
        Boolean disable_auto_fitting = false
        Boolean nodc = false
        Boolean no_ref = false
        Boolean just_ref = false
        Int threads = 1
        Int dcv = 1024
        Int offrate = 5
        Int ftabchars = 10
        Int localoffrate = 3
        Int localftabchars = 6
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(reference_fasta, "GiB") * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(reference_fasta, ".gz")}
        gunzip -c "~{reference_fasta}" > "$ref_fasta" \
            || ln -sf "~{reference_fasta}" "$ref_fasta"

        hisat2-build \
            ~{if force_large_index then "--large-index" else ""} \
            ~{if disable_auto_fitting then "--disable-auto-fitting" else ""} \
            -p ~{threads} \
            ~{if defined(bmax) then "--bmax \"~{bmax}\"" else ""} \
            ~{if defined(bmaxdivn) then "--bmaxdivn \"~{bmaxdivn}\"" else ""} \
            ~{if !nodc then "--dcv \"~{dcv}\"" else ""} \
            ~{if no_ref then "--no-ref" else ""} \
            ~{if just_ref then "--just-ref" else ""} \
            --offrate "~{offrate}" \
            --ftabchars "~{ftabchars}" \
            --localoffrate "~{localoffrate}" \
            --localftabchars "~{localftabchars}" \
            ~{if defined(snp) then "--snp \"~{snp}\"" else ""} \
            ~{if defined(haplotype) then "--haplotype \"~{haplotype}\"" else ""} \
            ~{if defined(splice_site) then "--ss \"~{splice_site}\"" else ""} \
            ~{if defined(exon) then "--exon \"~{exon}\"" else ""} \
            ~{if defined(repeat_ref) then "--repeat-ref \"~{repeat_ref}\"" else ""} \
            ~{if defined(repeat_info) then "--repeat-info \"~{repeat_info}\"" else ""} \
            ~{if defined(repeat_snp) then "--repeat-snp \"~{repeat_snp}\"" else ""} \
            ~{if defined(repeat_haplotype) then "--repeat-haplotype \"~{repeat_haplotype}\"" else ""} \
            ~{if defined(seed) then "--seed \"~{seed}\"" else ""} \
            "$ref_fasta" \
            "~{index_base_name}"

            tar -czf "~{index_base_name}.tar.gz" "~{index_base_name}"*
    >>>

    output {
        File reference_index = "~{index_base_name}.tar.gz"
    }

    requirements {
        cpu: threads
        memory: "16 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/hisat2:branch-minimap2-2.2.1-0"
    }
}
