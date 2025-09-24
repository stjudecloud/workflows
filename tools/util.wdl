## # Utilities

version 1.1

task download {
    meta {
        description: "Uses wget to download a file from a remote URL to the local filesystem"
        outputs: {
            downloaded_file: "File downloaded from provided URL"
        }
    }

    parameter_meta {
        url: "URL of the file to download"
        outfile_name: "Name of the output file"
        disk_size_gb: "Disk space to allocate for task, specified in GB"
        md5sum: "Optional md5sum to check against downloaded file. Recommended to use in order to catch corruption or an unintentional file swap."
    }

    input {
        String url
        String outfile_name
        Int disk_size_gb
        String? md5sum
    }

    command <<<
        set -euo pipefail

        wget "~{url}" -O "~{outfile_name}"

        if [ -n "~{md5sum}" ]; then
            echo "~{md5sum}  ~{outfile_name}" > "~{outfile_name}.md5"
            md5sum -c "~{outfile_name}.md5"
        fi
    >>>

    output {
        File downloaded_file = outfile_name
    }

    runtime {
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:3.0.0"
        maxRetries: 1
    }
}

task split_string {
    # Currently (v1.1) no way to do this using the WDL standard library.
    # Revisit task in future version updates, can hopefully be replaced.
    meta {
        description: "Split a string into an array of strings based on a delimiter"
        warning: "This implementation will result in a runtime error if the provided string has any embedded single quotes (`'`)!"
        outputs: {
            split_strings: "Split string as an array"
        }
    }

    parameter_meta {
        string: "String to split on occurences of `delimiter`"
        delimiter: {
            description: "Delimiter on which to split `input_string`",
        }
    }

    input {
        String string
        String delimiter
    }

    command <<<
        set -euo pipefail

        echo '~{sub(string, delimiter, "\n")}' > split_strings.txt
    >>>

    output {
        Array[String] split_strings = read_lines("split_strings.txt")
    }

    runtime {
        container: "ghcr.io/stjudecloud/util:3.0.0"
        maxRetries: 1
    }
}

task calc_feature_lengths {
    meta {
        description: "Calculate gene lengths from a GTF feature file using the non-overlapping exonic length algorithm"
        help: "The non-overlapping exonic length algorithm can be implemented as the sum of each base covered by at least one exon; where each base is given a value of 1 regardless of how many exons overlap it."
        outputs: {
            gene_lengths: "A two column headered TSV file with gene names in the first column and feature lengths (as integers) in the second column"
        }
    }

    parameter_meta {
        gtf: "GTF feature file"
        outfile_name: "Name of the gene lengths file"
        idattr: {
            description: "GTF attribute to be used as feature ID. The value of this attribute will be used as the first column in the output file.",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File gtf
        String outfile_name = basename(gtf, ".gtf.gz") + ".genelengths.txt"
        String idattr = "gene_name"
        Int modify_disk_size_gb = 0
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size_gb = ceil(gtf_size * 2) + 10 + modify_disk_size_gb

    command <<<
        python3 /scripts/util/calc_feature_lengths.py \
            --id_attr "~{idattr}" \
            "~{gtf}" \
            "~{outfile_name}"
    >>>

    output {
        File gene_lengths = outfile_name
    }

    runtime {
        memory: "16 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:3.0.0"
        maxRetries: 1
    }
}

task compression_integrity {
    meta {
        description: "Checks the compression integrity of a bgzipped file"
    }

    parameter_meta {
        bgzipped_file: "Input bgzipped file to check integrity of"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bgzipped_file
        Int modify_disk_size_gb = 0
    }

    Float file_size = size(bgzipped_file, "GiB")
    Int disk_size_gb = ceil(file_size) + 10 + modify_disk_size_gb

    command <<<
        bgzip -t "~{bgzipped_file}"
    >>>

    runtime {
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task add_to_bam_header {
    meta {
        description: "Adds another line of text to the bottom of a BAM header"
        outputs: {
            reheadered_bam: "The BAM after its header has been modified"
        }
    }

    parameter_meta {
        bam: "Input BAM format file which will have its header added to"
        additional_header: {
            description: "A string to add as a new line in the BAM header.",
            warning: "No format checking is done, so please ensure you do not invalidate your BAM with this task. Add only spec compliant entries to the header.",
        }
        prefix: "Prefix for the reheadered BAM. The extension `.bam` will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String additional_header
        String prefix = basename(bam, ".bam") + ".reheader"
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        samtools view -H "~{bam}" > header.sam
        echo "~{additional_header}" >> header.sam
        samtools reheader -P header.sam "~{bam}" > "~{outfile_name}"
    >>>

    output {
        File reheadered_bam = outfile_name
    }

    runtime {
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task unpack_tarball {
    meta {
        description: "Accepts a `.tar.gz` archive and converts it into a flat array of files. Any directory structure of the archive is ignored."
        outputs: {
            tarball_contents: "An array of files found in the input tarball"
        }
    }

    parameter_meta {
        tarball: "A `.tar.gz` archive to unpack into individual files"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File tarball
        Int modify_disk_size_gb = 0
    }

    Float tarball_size = size(tarball, "GiB")
    Int disk_size_gb = ceil(tarball_size * 8) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        mkdir unpacked_tarball
        tar -C unpacked_tarball -xzf "~{tarball}" --no-same-owner --no-same-permissions
        # pipe through sort because otherwise order is random (dependent on filesystem)
        find unpacked_tarball/ -type f | LC_ALL=C sort > file_list.txt
    >>>

    output {
        Array[File] tarball_contents = read_lines("file_list.txt")
    }

    runtime {
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:3.0.0"
        maxRetries: 1
    }
}

task make_coverage_regions_bed {
    meta {
        description: "Takes in a GTF file, converts it to BED, then filters it down to a 3 column BED file from all lines which match a given feature type"
        outputs: {
            bed: "3 column BED file corresponding to all records in the input GTF with a feature type matching `feature_type`",
        }
    }

    parameter_meta {
        gtf: "Gzipped GTF feature file from which to derive a coverage regions BED file"
        feature_type: {
            description: "Feature type to filter on. Only lines with this feature type will be included in the output BED file.",
            help: "`choices` below are the possible values from a GENCODE GTF file. If you are using a different GTF source, you may need to adjust this parameter.",
            choices: [
                "gene",
                "transcript",
                "exon",
                "CDS",
                "UTR",
                "start_codon",
                "stop_codon",
                "Selenocysteine",
            ],
        }
        outfile_name: "Name of the output BED file"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File gtf
        String feature_type
        String outfile_name = basename(gtf, "gtf.gz") + feature_type + ".bed"
        Int modify_disk_size_gb = 0
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size_gb = ceil(gtf_size * 1.2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        gunzip -c "~{gtf}" \
            | gtf2bed \
            | awk '$8 == "~{feature_type}" {print $1 "\t" $2 "\t" $3}' \
            > "~{outfile_name}"
    >>>

    output {
        File bed = outfile_name
    }

    runtime {
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0"
        maxRetries: 1
    }
}

task global_phred_scores {
    meta {
        description: "Calculates statistics about PHRED scores of the input BAM"
        outputs: {
            phred_scores: "Headered TSV file containing PHRED score statistics"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to calculate PHRED score statistics for"
        prefix: "Prefix for the output TSV file. The extension `.global_PHRED_scores.tsv` will be added."
        fast_mode: "Enable fast mode (true) or calculate statistics for *_every_* base in the BAM (false)?"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean fast_mode = true
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".global_PHRED_scores.tsv"

    command <<<
        python3 /scripts/util/calc_global_phred_scores.py \
            ~{if fast_mode then "--fast_mode" else ""} \
            "~{bam}" \
            "~{prefix}"
    >>>

    output {
        File phred_scores = "~{outfile_name}"
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:3.0.0"
        maxRetries: 1
    }
}

task check_fastq_and_rg_concordance {
    meta {
        description: "Validates FASTQs and read group records are concordant"
        help: "Each read1 FASTQ must correspond to exactly one read group record. This correspondance is encoded in two ways, both of which must match. 1) the FASTQ and its read group share the same index of their respective lists and 2) the `ID` field value must be contained somewhere within the FASTQ file basename. If `read_two_names` is non-empty, the same checks are performed on each of these names as well (i.e. all 3 of the read1 FASTQ, read2 FASTQ, and read group ID must match and be in the same position of their list). Additionally, the `ID` field must be the first field of the read group record, and each `ID` value must be unique."
        warning: "This task does not do any validation outside of what is described. i.e. only the first field of read group records are checked, and any malformed records beyond that field will go undetected."
    }

    parameter_meta {
        read_one_names: {
            description: "Filenames of every read1 FASTQ to validate.",
            help: "Either basenames or full paths are acceptable. Must be non-empty.",
        }
        read_groups: {
            description: "Read group records that correspond to the FASTQs being validated.",
            help: "Read group records may be optionally prefixed with `@RG` and may use either escaped tabs (`\t`) or spaces as delimiters.",
        }
        read_two_names: {
            description: "Filenames of every read2 FASTQ to validate.",
            help: "Either basenames or full paths are acceptable. May be empty or the exact same length as `read_one_names` and `read_groups`.",
        }
    }

    input {
        Array[String] read_one_names
        Array[String] read_groups
        Array[String]? read_two_names
    }

    Array[String] read_twos = select_first([read_two_names, []])

    command <<<
        python3 /scripts/util/check_FQs_and_RGs.py \
            --read-one-fastqs "~{sep(",", read_one_names)}" \
            ~{(
                if length(read_twos) > 0
                then "--read-two-fastqs \"" + sep(",", squote(read_twos)) + "\""
                else ""
            )} \
            --read-groups "~{sep(",", read_groups)}"
    >>>

    runtime {
        container: "ghcr.io/stjudecloud/util:3.0.0"
        maxRetries: 1
    }
}

task split_fastq {
    meta {
        description: "Splits a FASTQ into multiple files based on the number of reads per file"
        outputs: {
            fastqs: "Array of FASTQ files, each containing a subset of the input FASTQ"
        }
    }

    parameter_meta {
        fastq: {
            description: "Gzipped FASTQ file to split",
            stream: true,
        }
        prefix: {
            description: "Prefix for the FASTQ files. The extension `.fastq.gz` (preceded by a split index) will be added.",
            group: "Common",
        }
        reads_per_file: "Number of reads to include in each output FASTQ file"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File fastq
        String prefix = sub(
            basename(fastq),
            "(fastq|fq)\\.gz$",
            ""
        )
        Int reads_per_file = 10000000
        Int modify_disk_size_gb = 0
        Int ncpu = 2
    }

    Float fastq_size = size(fastq, "GiB")
    Int disk_size_gb = ceil(fastq_size * 5) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        (( lines = ~{reads_per_file} * 4 ))
        zcat "~{fastq}" | split -l $lines -d -a 6 - "~{prefix}"

        for file in "~{prefix}"*; do
            mv "$file" "${file}.fastq"
            echo "gzip ${file}.fastq" >> cmds
        done

        parallel --jobs ~{ncpu} < cmds
    >>>

    output {
        Array[File] fastqs = glob("~{prefix}*")
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:3.0.0"
        maxRetries: 1
    }
}
