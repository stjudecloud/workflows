version 1.1

import "../../data_structures/read_group.wdl"
import "./fastq-to-ubam-core.wdl" as ubam_core

workflow fastq_to_ubam {
    meta {
        description: "Generate unaligned BAM file from FASTQ files"
        outputs: {
            unaligned_bam: "Unaligned BAM file",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        read_one_fastqs_gz: "Array of gzipped FASTQ files with 1st reads in pair"
        read_two_fastqs_gz: "Array of gzipped FASTQ files with 2nd reads in pair"
        read_groups: {
            description: "An Array of structs defining read groups to include in the harmonized BAM. Must correspond to input FASTQs. Each read group ID must be contained in the basename of a FASTQ file or pair of FASTQ files if Paired-End. This requirement means the length of `read_groups` must equal the length of `read_one_fastqs_gz` and the length of `read_two_fastqs_gz` if non-zero. Only the `ID` field is required, and it must be unique for each read group defined. See top of file for help formatting your input JSON.",  # Does not currently handle unknown RG case
            help: "See ../../data_structures/read_group.wdl for additional documentation.",
            external_help: "https://samtools.github.io/hts-specs/SAMv1.pdf",
        }
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        validate_input: "Ensure input BAM is well-formed before beginning harmonization?"
        use_all_cores: "Use all cores? Recommended for cloud environments."
    }

    input {
        Array[File] read_one_fastqs_gz
        Array[File] read_two_fastqs_gz
        Array[ReadGroup] read_groups
        String prefix = basename(read_one_fastqs_gz[0], ".fastq.gz")
        Boolean validate_input = true
        Boolean use_all_cores = false
    }

    call parse_input { input:
        read_one_fastqs = length(read_one_fastqs_gz),
        read_two_fastqs = length(read_two_fastqs_gz),
        read_groups,
    }

    scatter (read_group in read_groups) {
        call read_group.validate_read_group as validate_readgroups { input:
            read_group,
            required_fields = ["ID", "SM", "LB", "CN", "DT"]
        }
    }

    call ubam_core.fastq_to_ubam_core after parse_input after validate_readgroups { input:
        read_one_fastqs_gz,
        read_two_fastqs_gz,
        read_groups,
        prefix,
        validate_input,
        use_all_cores,
    }

    output {
        File unaligned_bam = fastq_to_ubam_core.unaligned_bam
    }
}

task parse_input {
    meta {
        description: "Parses and validates the `fastq-to-ubam` workflow's provided inputs"
        outputs: {
            check: "Dummy output to indicate success and to enable call-caching"
        }
    }

    parameter_meta {
        read_groups: {
            description: "An Array of structs defining read groups to include in the harmonized BAM. Must correspond to input FASTQs. Each read group ID must be contained in the basename of a FASTQ file or pair of FASTQ files if Paired-End. This requirement means the length of `read_groups` must equal the length of `read_one_fastqs_gz` and the length of `read_two_fastqs_gz` if non-zero. Only the `ID` field is required, and it must be unique for each read group defined. See top of file for help formatting your input JSON.",  # Does not currently handle unknown RG case
            help: "See ../../data_structures/read_group.wdl for additional documentation.",
            external_help: "https://samtools.github.io/hts-specs/SAMv1.pdf",
        }
        read_one_fastqs: "Number of gzipped FASTQ files with 1st reads in pair"
        read_two_fastqs: "Number of gzipped FASTQ files with 2nd reads in pair"
    }

    input {
        Array[ReadGroup] read_groups
        Int read_one_fastqs
        Int read_two_fastqs
    }

    #@ except: LineWidth
    command <<<
        if [ ~{read_one_fastqs} -ne ~{read_two_fastqs} ]
        then
            >&2 echo "Number of entries in read_one_fastqs_gz and read_two_fastqs_gz must be equal"
            exit 1
        fi
        if [ ~{read_one_fastqs} -ne ~{length(read_groups)} ]
        then
            >&2 echo "Number of entries read_one_fastqs_gz and read_groups must be equal"
            exit 1
        fi
        if [ ~{read_two_fastqs} -ne ~{length(read_groups)} ]
        then
            >&2 echo "Number of entries read_two_fastqs_gz and read_groups must be equal"
            exit 1
        fi
    >>>

    output {
        String check = "passed"
    }

    runtime {
        memory: "4 GB"
        disks: "10 GB"
        container: "ghcr.io/stjudecloud/util:2.1.0"
        maxRetries: 1
    }
}
