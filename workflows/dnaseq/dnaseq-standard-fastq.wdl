## **WARNING:** this workflow is experimental! Use at your own risk!

version 1.1

import "../../data_structures/read_group.wdl"
import "../../tools/fq.wdl"
import "./dnaseq-core.wdl" as dnaseq_core_wf

workflow dnaseq_standard_fastq_experimental {
    meta {
        description: "Aligns DNA reads using bwa"
        outputs: {
            harmonized_bam: "Harmonized DNA-Seq BAM, aligned with bwa",
            harmonized_bam_index: "Index for the harmonized DNA-Seq BAM file"
        }
        allowNestedInputs: true
    }
    parameter_meta {
        read_one_fastqs_gz: "Input gzipped FASTQ format file(s) with 1st read in pair to align"
        read_two_fastqs_gz: "Input gzipped FASTQ format file(s) with 2nd read in pair to align"
        bwa_db: "Gzipped tar archive of the bwa reference files. Files should be at the root of the archive."
        reads_per_file: "Number of reads per FASTQ file to output."
        read_groups: {
            description: "An Array of structs defining read groups to include in the harmonized BAM. Must correspond to input FASTQs. Each read group ID must be contained in the basename of a FASTQ file or pair of FASTQ files if Paired-End. This requirement means the length of `read_groups` must equal the length of `read_one_fastqs_gz` and the length of `read_two_fastqs_gz` if non-zero. Only the `ID` field is required, and it must be unique for each read group defined. See data_structures/read_group.wdl for help formatting your input JSON.",
            external_help: "https://samtools.github.io/hts-specs/SAMv1.pdf",
        }
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        aligner: {
            description: "BWA aligner to use",
            choices: ["mem", "aln"]
        }
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
        subsample_n_reads: "Only process a random sampling of `n` reads. Any `n`<=`0` for processing entire input."
    }
    input {
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        File bwa_db
        Int reads_per_file = 10000000
        Array[ReadGroup] read_groups
        String prefix
        String aligner = "mem"
        Boolean validate_input = true
        Boolean use_all_cores = false
        Int subsample_n_reads = -1
    }

    call dnaseq_core_wf.parse_input { input:
        aligner
    }

    if (validate_input){
        scatter (reads in zip(read_one_fastqs_gz, read_two_fastqs_gz)) {
            call fq.fqlint { input:
                read_one_fastq=reads.left,
                read_two_fastq=reads.right,
            }
        }
    }

    if (subsample_n_reads > 0) {
        Int reads_per_pair = ceil(subsample_n_reads / length(read_one_fastqs_gz))
        scatter (reads in zip(read_one_fastqs_gz, read_two_fastqs_gz)) {
            call fq.subsample { input:
                read_one_fastq=reads.left,
                read_two_fastq=reads.right,
                record_count=reads_per_pair,
            }
        }
    }
    Array[File] selected_read_one_fastqs = select_first([
        subsample.subsampled_read1,
        read_one_fastqs_gz
    ])
    Array[File] selected_read_two_fastqs = select_all(
        select_first([
            subsample.subsampled_read2,
            read_two_fastqs_gz
        ])
    )

    call dnaseq_core_wf.dnaseq_core_experimental { input:
        read_one_fastqs_gz=selected_read_one_fastqs,
        read_two_fastqs_gz=selected_read_two_fastqs,
        bwa_db=bwa_db,
        reads_per_file=reads_per_file,
        read_groups=read_groups,
        prefix=prefix,
        aligner=aligner,
        use_all_cores=use_all_cores
    }

    output {
        File harmonized_bam = dnaseq_core_experimental.harmonized_bam
        File harmonized_bam_index = dnaseq_core_experimental.harmonized_bam_index
    }
}
