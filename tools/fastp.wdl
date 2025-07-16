version 1.1

task fastp {
    meta {
        description: "Runs the Fastp tool for FASTQ quality control and trimming"
        outputs: {
            fastp_report: "A gzipped tar archive of all Fastp output files"
        }
    }

    parameter_meta {
        read_one_fastq: "Input FASTQ with read one. Can be gzipped or uncompressed."
        read_two_fastq: "Input FASTQ with read two. Can be gzipped or uncompressed."
        prefix: "Prefix for the Fastp report files. The extensions `.fastp.html`, `.fastp.json`, and `.fastp.trimmed.fastq.gz` will be added."
        disable_adapter_trimming: "Disable adapter trimming"
        deduplicate: "Enable deduplication to drop the duplicated reads/pairs"
        phred64: "Input uses phred64 encoding. It will be converted to phred33 encoding in the output files."
        trim_front_r1: "Number of bases to trim from the front of read one"
        trim_tail_r1: "Number of bases to trim from the tail of read one"
        trim_front_r2: "Number of bases to trim from the front of read two"
        trim_tail_r2: "Number of bases to trim from the tail of read two"
        max_length_r1: "Maximum length of read one. Reads longer than this will be trimmed from the tail."
        max_length_r2: "Maximum length of read two. Reads longer than this will be trimmed from the tail."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File read_one_fastq
        File? read_two_fastq
        String prefix = sub(
            basename(read_one_fastq),
            "(([_.][rR](?:ead)?[12])((?:[_.-][^_.-]*?)*?))?\\.(fastq|fq)(\\.gz)?$",
            ""  # Once replacing with capturing groups is supported, replace with group 3
        )
        Boolean disable_adapter_trimming = false
        Boolean deduplicate = false
        Boolean phred64 = false
        Int trim_front_r1 = 0
        Int trim_tail_r1 = 0
        Int trim_front_r2 = 0
        Int trim_tail_r2 = 0
        Int max_length_r1 = 0
        Int max_length_r2 = 0
        Int modify_disk_size_gb = 0
        Int ncpu = 3
    }

    Float input_size = size(read_one_fastq, "GiB")
        + (if defined(read_two_fastq) then size(read_two_fastq, "GiB") else 0)
    Int disk_size_gb = ceil(input_size) + 10 + modify_disk_size_gb

    String out_tar_gz = prefix + ".tar.gz"

    command <<< 
        set -euo pipefail

        # set ENV variables for `fastp`
        export LC_ALL=C.UTF-8
        export LANG=C.UTF-8

        fastp \
            -i "~{read_one_fastq}" \
            ~{"-I '" + read_two_fastq + "'"} \
            -o "~{prefix}.R1.fastq.gz" \
            ~{if defined(read_two_fastq) then "-O '" + prefix + ".R2.fastq.gz'" else ""} \
            ~{if phred64 then "--phred64" else ""} \
            ~{if disable_adapter_trimming then "--disable_adapter_trimming" else ""} \
            --trim_front1 ~{trim_front_r1} \
            --trim_tail1 ~{trim_tail_r1} \
            --trim_front2 ~{trim_front_r2} \
            --trim_tail2 ~{trim_tail_r2} \
            --max_len1 ~{max_length_r1} \
            --max_len2 ~{max_length_r2} \
            -R "~{prefix} report" \
            --thread ~{ncpu} \
            ~{if deduplicate then "--dedup" else ""} \

            -h "~{prefix}.fastp.html" \
            -j "~{prefix}.fastp.json"

        tar -czf ~{out_tar_gz} \
            ~{prefix}.R1.fastq.gz \
            ~{prefix}.R2.fastq.gz \
            ~{prefix}.fastp.html \
            ~{prefix}.fastp.json
    >>>

    output {
        File fastp_report = out_tar_gz
    }

    runtime {
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/fastp:1.0.1--heae3180_0"
        maxRetries: 1
    }
}
