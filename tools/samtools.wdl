## [Homepage](http://samtools.sourceforge.net/)
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
version 1.1

task quickcheck {
    meta {
        description: "Runs Samtools quickcheck on the input BAM file. This checks that the BAM file appears to be intact, e.g. header exists, at least one sequence is present, and the end-of-file marker exists."
        outputs: {
            checked_bam: "The unmodfied input BAM after it has been successfully quickchecked"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to quickcheck"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        samtools quickcheck ~{bam}
    >>>

    output {
        File checked_bam = bam
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task split {
    meta {
        description: "Runs Samtools split on the input BAM file. This splits the BAM by read group into one or more output files. It optionally errors if there are reads present that do not belong to a read group."
    }

    parameter_meta {
        bam: "Input BAM format file to split"
        prefix: "Prefix for the split BAM files. The extensions will contain read group IDs, and will end in `.bam`."
        reject_unaccounted: {
            description: "If true, error if there are reads present that do not have read group information.",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean reject_unaccounted = true
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools split \
            --threads "$n_cores" \
            -u ~{prefix}.unaccounted_reads.bam \
            -f '~{prefix}.%!.bam' \
            ~{bam}

        samtools head \
            --threads "$n_cores" \
            -h 0 \
            -n 1 \
            ~{prefix}.unaccounted_reads.bam \
            > first_unaccounted_read.bam

        if ~{reject_unaccounted} && [ -s first_unaccounted_read.bam ]; then
            exit 42
        else
            rm ~{prefix}.unaccounted_reads.bam
        fi
        rm first_unaccounted_read.bam
    >>>

    output {
       Array[File] split_bams = glob("*.bam")
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task flagstat {
    meta {
        description: "Produces a `samtools flagstat` report containing statistics about the alignments based on the bit flags set in the BAM"
        outputs: {
            flagstat_report: "`samtools flagstat` STDOUT redirected to a file"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to generate flagstat for"
        outfile_name: "Name for the flagstat report file"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".flagstat.txt"
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
        # TODO should we support alt formats? Will they break MultiQC?
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools flagstat --threads "$n_cores" ~{bam} > ~{outfile_name}
    >>>

    output {
       File flagstat_report = outfile_name
    }

    runtime {
        memory: "5 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task index {
    meta {
        description: "Creates a `.bai` BAM index for the input BAM"
        outputs: {
            bam_index: "A `.bai` BAM index associated with the input BAM. Filename will be `basename(bam) + '.bai'`."
        }
    }

    parameter_meta {
        bam: "Input BAM format file to index"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 1.2) + 10 + modify_disk_size_gb

    String outfile_name = basename(bam) + ".bai"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        # For some reason, index doesn't support '--threads',
        # so we use '-@' here
        samtools index -@ "$n_cores" ~{bam} ~{outfile_name}
    >>>

    output {
        File bam_index = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task subsample {
    meta {
        description: "Randomly subsamples the input BAM. Sampling is probabalistic and will be approximate to `desired_reads`. Read count will not be exact. A `sampled_bam` will not be produced if the input BAM read count is less than or equal to `desired_reads`."
        outputs: {
            orig_read_count: "A TSV report containing the original read count before subsampling",
            sampled_bam: "The subsampled input BAM."
        }
    }

    parameter_meta {
        bam: "Input BAM format file to subsample"
        desired_reads: "How many reads should be in the ouput BAM? Output BAM read count will be approximate to this value."
        prefix: "Prefix for the BAM file. The extension `.subsampled.bam` will be added."
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        Int desired_reads
        String prefix = basename(bam, ".bam")
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    String suffixed = prefix + ".subsampled"

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        if [[ ~{desired_reads} -le 0 ]]; then
            echo "'desired_reads' must be >0!" > /dev/stderr
            exit 1
        fi

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        read_count="$(samtools view --threads "$n_cores" -c ~{bam})"

        if [[ "$read_count" -gt "~{desired_reads}" ]]; then
            # the BAM has at least ~{desired_reads} reads, meaning we should
            # subsample it.
            frac=$( \
                awk -v desired_reads=~{desired_reads} \
                    -v read_count="$read_count" \
                        'BEGIN{
                            printf "%1.8f",
                            ( desired_reads / read_count )
                        }' \
                )
            samtools view --threads "$n_cores" -hb -s "$frac" ~{bam} \
                > ~{suffixed}.bam

            {
                echo -e "sample\toriginal read count"
                echo -e "~{prefix}\t$read_count"
            } > ~{suffixed}.orig_read_count.tsv
        else
            # the BAM has less than ~{desired_reads} reads, meaning we should
            # just use it directly without subsampling.

            # Do not use the '.subsampled' suffixed name
            # if not subsampled. Use ~{prefix} instead.
            {
                echo -e "sample\toriginal read count"
                echo -e "~{prefix}\t-"
            } > ~{prefix}.orig_read_count.tsv
        fi
    >>>

    output {
        File orig_read_count = glob("*.orig_read_count.tsv")[0]
        File? sampled_bam = suffixed + ".bam"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task merge {
    meta {
        description: "Merges multiple sorted BAMs into a single BAM"
        outputs: {
            merged_bam: "The BAM resulting from merging all the input BAMs"
        }
    }

    parameter_meta {
        bams: "An array of BAMs to merge into one combined BAM"
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        new_header: "Use the lines of FILE as `@` headers to be copied to the merged BAM, replacing any header lines that would otherwise be copied from the first BAM file in the list. (File may actually be in SAM format, though any alignment records it may contain are ignored.)"
        region: "Merge files in the specified region (Format: `chr:start-end`)"
        attach_rg: {
            description: "Attach an RG tag to each alignment. The tag value is inferred from file names.",
            common: true
        }
        name_sorted: {
            description: "Are _all_ input BAMs `queryname` sorted (true)? Or are _all_ input BAMs `coordinate` sorted (false)?",
            common: true
        }
        combine_rg: {
            description: "When several input files contain @RG headers with the same ID, emit only one of them (namely, the header line from the first file we find that ID in) to the merged output file. Combining these similar headers is usually the right thing to do when the files being merged originated from the same file. Without `-c`, all @RG headers appear in the output file, with random suffixes added to their IDs where necessary to differentiate them.",
            common: true
        }
        combine_pg: {
            description: "Similarly to `combine_rg`: for each @PG ID in the set of files to merge, use the @PG line of the first file we find that ID in rather than adding a suffix to differentiate similar IDs.",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] bams
        String prefix
        File? new_header
        String region = ""
        Boolean attach_rg = true
        Boolean name_sorted = false
        Boolean combine_rg = true
        Boolean combine_pg = true
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bams_size = size(bams, "GiB")
    Float header_size = size(new_header, "GiB")
    Int disk_size_gb = ceil(bams_size * 2 + header_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools merge \
            --threads "$n_cores" \
            ~{if defined(new_header) then "-h " + new_header else ""} \
            ~{if name_sorted then "-n" else ""} \
            ~{if (region != "") then "-R " + region else ""} \
            ~{if attach_rg then "-r" else ""} \
            ~{if combine_rg then "-c" else ""} \
            ~{if combine_pg then "-p" else ""} \
            ~{prefix}.bam \
            ~{sep(" ", bams)}
    >>>

    output {
        File merged_bam = prefix + ".bam"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task addreplacerg {
    meta {
        description: "Adds or replaces read group tags"
        outputs: {
            tagged_bam: "The transformed input BAM after read group modifications have been applied"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to add read group information"
        read_group_id: "Allows you to specify the read group ID of an existing @RG line and applies it to the reads specified by the `orphan_only` option"
        read_group_line: {
            description: "Allows you to specify a read group line to append to (or replace in) the header and applies it to the reads specified by the `orphan_only` option. Each String in the Array should correspond to one field of the read group line. Tab literals will be inserted between each entry in the final BAM. Only **one** read group line can be supplied per invocation of this task.",
            common: true
        }
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        orphan_only: {
            description: "Only add RG tags to orphans (true)? Or _also_ overwrite all existing RG tags (including any in the header) (false)?",
            common: true
        }
        overwrite_header_record: {
            description: "Overwrite an existing @RG line, if a new one with the same ID value is provided?",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String? read_group_id
        Array[String] read_group_line = []
        String prefix = basename(bam, ".bam") + ".addreplacerg"
        Boolean orphan_only = true
        Boolean overwrite_header_record = false
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools addreplacerg \
            --threads "$n_cores" \
            ~{sep(" ", prefix("-r ", squote(read_group_line)))} \
            ~{if defined(read_group_id) then "-R " + read_group_id else ""} \
            -m ~{if orphan_only then "orphan_only" else "overwrite_all"} \
            ~{if overwrite_header_record then "-w" else ""} \
            -o ~{outfile_name} \
            ~{bam}
    >>>

    output {
        File tagged_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task collate {
    meta {
        description: "Runs `samtools collate` on the input BAM file. Shuffles and groups reads together by their names."
        outputs: {
            collated_bam: "A collated BAM (reads sharing a name next to each other, no other guarantee of sort order)"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to collate"
        prefix: "Prefix for the collated BAM file. The extension `.bam` will be added."
        fast_mode: {
            description: "Use fast mode (output primary alignments only)?",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".collated"
        Boolean fast_mode = true
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = ceil(bam_size * 0.2) + 4 + modify_memory_gb
    Int disk_size_gb = ceil(bam_size * 4) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools collate \
            --threads "$n_cores" \
            ~{if fast_mode then "-f" else ""} \
            -o ~{outfile_name} \
            ~{bam}
    >>>

    output {
        File collated_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task bam_to_fastq {
    meta {
        description: "Runs `samtools fastq` on the input BAM file. Converts the BAM into FASTQ files. If `paired_end = false`, then _all_ reads in the BAM will be output to a single FASTQ file. Use filtering arguments to remove any unwanted reads. Assumes either a name sorted or collated BAM. For splitting a position sorted BAM see `collate_to_fastq`."
        outputs: {
            read_one_fastq_gz: "Gzipped FASTQ file with 1st reads in pair",
            read_two_fastq_gz: "Gzipped FASTQ file with 2nd reads in pair",
            singleton_reads_fastq_gz: "A gzipped FASTQ containing singleton reads",
            interleaved_reads_fastq_gz: "An interleaved gzipped paired-end FASTQ",
            single_end_reads_fastq_gz: "A gzipped FASTQ containing all reads"
        }
    }

    parameter_meta {
        bam: "Input name sorted or collated BAM format file to convert into FASTQ(s)"
        prefix: "Prefix for output FASTQ(s). Extensions `[,.R1,.R2,.singleton].fastq.gz` will be added depending on other options."
        f: "Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        F: "Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/). This defaults to 0x900 representing filtering of secondary and supplementary alignments."
        # rf: "Only output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/), in octal by beginning with `0` (i.e. /^0[0-7]+/)."  # introduced in v1.18 no quay.io image yet
        G: "Only EXCLUDE reads with all of the bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        paired_end: {
            description: "Is the data Paired-End? If `paired_end = false`, then _all_ reads in the BAM will be output to a single FASTQ file. Use filtering arguments to remove any unwanted reads.",
            common: true
        }
        append_read_number: {
            description: "Append /1 and /2 suffixes to read names",
            common: true
        }
        interleaved: {
            description: "Create an interleaved FASTQ file from Paired-End data? Ignored if `paired_end = false`.",
            common: true
        }
        output_singletons: "Output singleton reads as their own FASTQ?"
        fail_on_unexpected_reads: {
            description: "Should the task fail if reads with an unexpected `first`/`last` bit setting are discovered?",
            help: "The definition of 'unexpected' depends on whether the values of `paired_end` and `output_singletons` are true or false. In any case, reads that have neither or both `first` and `last` bits set are considered unexpected. If `paired_end` is `true` and `output_singletons` is `false`, singleton reads are considered unexpected. A singleton read is a read with either the `first` or the `last` bit set and that possesses a _unique_ QNAME; i.e. it is a read without a pair when all reads are expected to be paired. But if `output_singletons` is `true`, these singleton reads will be output as their own FASTQ instead of causing the task to fail. If `fail_on_unexpected_reads` is `false`, then all the above cases will be ignored. Any 'unexpected' reads will be silently discarded. If `paired_end` is `false`, no reads are considered unexpected, and _every_ read will be present in the resulting FASTQ regardless of bit settings.",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        String f = "0"
        String F = "0x900"
        # String rf = "0"  # introduced in v1.18 no quay.io image yet
        String G = "0"
        Boolean paired_end = true
        Boolean append_read_number = true
        Boolean interleaved = false
        Boolean output_singletons = false
        Boolean fail_on_unexpected_reads = false
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools fastq \
            --threads "$n_cores" \
            -f ~{f} \
            -F ~{F} \
            -G ~{G} \
            ~{if append_read_number
                then "-N"
                else "-n"
            } \
            -1 ~{
                if paired_end then (
                    if interleaved then prefix + ".fastq.gz" else prefix + "_R1.fastq.gz"
                )
                else prefix + ".fastq.gz"
            } \
            -2 ~{
                if paired_end then (
                    if interleaved then prefix + ".fastq.gz" else prefix + ".R2.fastq.gz"
                )
                else prefix + ".fastq.gz"
            } \
            ~{
                if paired_end then (
                    if output_singletons
                    then "-s " + prefix+".singleton.fastq.gz"
                    else "-s junk.singleton.fastq.gz"
                )
                else ""
            } \
            -0 ~{
                if paired_end
                then "junk.unknown_bit_setting.fastq.gz"
                else prefix + ".fastq.gz"
            } \
            ~{bam}

        if ~{fail_on_unexpected_reads} \
            && find . -name 'junk.*.fastq.gz' ! -empty | grep -q .
        then
            >&2 echo "Discovered unexpected reads in:"
            find . -name 'junk.*.fastq.gz' ! -empty >&2
            exit 42
        fi
    >>>

    output {
        File? read_one_fastq_gz = "~{prefix}.R1.fastq.gz"
        File? read_two_fastq_gz = "~{prefix}.R2.fastq.gz"
        File? singleton_reads_fastq_gz = "~{prefix}.singleton.fastq.gz"
        File? interleaved_reads_fastq_gz = "~{prefix}.fastq.gz"
        File? single_end_reads_fastq_gz = "~{prefix}.fastq.gz"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task collate_to_fastq {
    meta {
        description: "Runs `samtools collate` on the input BAM file then converts it into FASTQ(s) using `samtools fastq`"
        outputs: {
            collated_bam: "A collated BAM (reads sharing a name next to each other, no other guarantee of sort order)",
            read_one_fastq_gz: "Gzipped FASTQ file with 1st reads in pair",
            read_two_fastq_gz: "Gzipped FASTQ file with 2nd reads in pair",
            singleton_reads_fastq_gz: "Gzipped FASTQ containing singleton reads",
            interleaved_reads_fastq_gz: "Interleaved gzipped paired-end FASTQ",
            single_end_reads_fastq_gz: "A gzipped FASTQ containing all reads"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to collate and convert to FASTQ(s)"
        prefix: "Prefix for the collated BAM and FASTQ files. The extensions `.collated.bam` and `[,.R1,.R2,.singleton].fastq.gz` will be added."
        f: "Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        F: "Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/). This defaults to 0x900 representing filtering of secondary and supplementary alignments."
        # rf: "Only output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/), in octal by beginning with `0` (i.e. /^0[0-7]+/)."  # introduced in v1.18 no quay.io image yet
        G: "Only EXCLUDE reads with all of the bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        fast_mode: {
            description: "Fast mode for `samtools collate` (primary alignments only)",
            common: true
        }
        store_collated_bam: {
            description: "Save the collated BAM to disk and output it (true)?",
            common: true
        }
        paired_end: {
            description: "Is the data Paired-End? If `paired_end = false`, then _all_ reads in the BAM will be output to a single FASTQ file. Use filtering arguments to remove any unwanted reads.",
            common: true
        }
        append_read_number: {
            description: "Append /1 and /2 suffixes to read names",
            common: true
        }
        interleaved: {
            description: "Create an interleaved FASTQ file from Paired-End data? Ignored if `paired_end = false`.",
            common: true
        }
        output_singletons: "Output singleton reads as their own FASTQ?"
        fail_on_unexpected_reads: {
            description: "Should the task fail if reads with an unexpected `first`/`last` bit setting are discovered?",
            help: "The definition of 'unexpected' depends on whether the values of `paired_end` and `output_singletons` are true or false. In any case, reads that have neither or both `first` and `last` bits set are considered unexpected. If `paired_end` is `true` and `output_singletons` is `false`, singleton reads are considered unexpected. A singleton read is a read with either the `first` or the `last` bit set and that possesses a _unique_ QNAME; i.e. it is a read without a pair when all reads are expected to be paired. But if `output_singletons` is `true`, these singleton reads will be output as their own FASTQ instead of causing the task to fail. If `fail_on_unexpected_reads` is `false`, then all the above cases will be ignored. Any 'unexpected' reads will be silently discarded. If `paired_end` is `false`, no reads are considered unexpected, and _every_ read will be present in the resulting FASTQ regardless of bit settings.",
            common: true
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        String f = "0"
        String F = "0x900"
        # String rf = "0"  # introduced in v1.18 no quay.io image yet
        String G = "0"
        Boolean fast_mode = true
        Boolean store_collated_bam = false
        Boolean paired_end = true
        Boolean append_read_number = true
        Boolean interleaved = false
        Boolean output_singletons = false
        Boolean fail_on_unexpected_reads = false
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = ceil(bam_size * 0.4) + 4 + modify_memory_gb
    Int disk_size_gb = ceil(bam_size * 5) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        # Use the `-u` flag to skip compression (and decompression)
        # if not storing the output
        samtools collate \
            ~{if store_collated_bam then "" else "-u"} \
            --threads "$n_cores" \
            ~{if fast_mode then "-f" else ""} \
            -O \
            ~{bam} \
            | tee ~{if store_collated_bam then prefix + ".collated.bam" else ""} \
            | samtools fastq \
                --threads "$n_cores" \
                -f ~{f} \
                -F ~{F} \
                -G ~{G} \
                ~{if append_read_number
                    then "-N"
                    else "-n"
                } \
                -1 ~{
                    if paired_end then (
                        if interleaved then prefix + ".fastq.gz" else prefix + ".R1.fastq.gz"
                    )
                    else prefix + ".fastq.gz"
                } \
                -2 ~{
                    if paired_end then (
                        if interleaved then prefix + ".fastq.gz" else prefix + ".R2.fastq.gz"
                    )
                    else prefix + ".fastq.gz"
                } \
                ~{
                    if paired_end then (
                        if output_singletons
                        then "-s " + prefix+".singleton.fastq.gz"
                        else "-s junk.singleton.fastq.gz"
                    )
                    else ""
                } \
                -0 ~{
                    if paired_end
                    then "junk.unknown_bit_setting.fastq.gz"
                    else prefix + ".fastq.gz"
                }

        if ~{fail_on_unexpected_reads} \
            && find . -name 'junk.*.fastq.gz' ! -empty | grep -q .
        then
            >&2 echo "Discovered unexpected reads in:"
            find . -name 'junk.*.fastq.gz' ! -empty >&2
            exit 42
        fi
    >>>

    output {
        File? collated_bam = "~{prefix}.collated.bam"
        File? read_one_fastq_gz = "~{prefix}.R1.fastq.gz"
        File? read_two_fastq_gz = "~{prefix}.R2.fastq.gz"
        File? singleton_reads_fastq_gz = "~{prefix}.singleton.fastq.gz"
        File? interleaved_reads_fastq_gz = "~{prefix}.fastq.gz"
        File? single_end_reads_fastq_gz = "~{prefix}.fastq.gz"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}

task fixmate {
    meta {
        description: "Runs `samtools fixmate` on the input BAM file. This fills in mate coordinates and insert size fields."
    }

    parameter_meta {
        bam: "Input BAM format file to add mate information. Must be name-sorted or name-collated."
        prefix: "Prefix for the output file. The extension specified with the `extension` parameter will be added."
        extension: {
            description: "File format extension to use for output file.",
            choices: [
                ".bam",
                ".cram"
            ],
            common: true
        }
        add_cigar: "Add template cigar ct tag"
        add_mate_score: "Add mate score tags. These are used by markdup to select the best reads to keep."
        disable_proper_pair_check: "Disable proper pair check [ensure one forward and one reverse read in each pair]"
        remove_unaligned_and_secondary: "Remove unmapped and secondary reads"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            common: true
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".fixmate"
        String extension = ".bam"
        Boolean add_cigar = true
        Boolean add_mate_score = true
        Boolean disable_proper_pair_check = false
        Boolean remove_unaligned_and_secondary = false
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools fixmate \
            --threads "$n_cores" \
            ~{if remove_unaligned_and_secondary then "-r" else ""} \
            ~{if disable_proper_pair_check then "-p" else ""} \
            ~{if add_cigar then "-c" else ""} \
            ~{if add_mate_score then "-m" else ""} \
            ~{bam} ~{prefix}~{extension}
    >>>

    output {
        File fixmate_bam = "~{prefix}~{extension}"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'docker://quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: 1
    }
}
