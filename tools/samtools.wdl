## [Homepage](http://samtools.sourceforge.net/)

version 1.1

import "../data_structures/flag_filter.wdl"

task quickcheck {
    meta {
        description: "Runs Samtools quickcheck on the input BAM file. This checks that the BAM file appears to be intact, e.g. header exists and the end-of-file marker exists."
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
        samtools quickcheck "~{bam}"
    >>>

    runtime {
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task split {
    meta {
        description: "Runs Samtools split on the input BAM file. This splits the BAM by read group into one or more output files. It optionally errors if there are reads present that do not belong to a read group."
        outputs: {
            split_bams: "The split BAM files. The extensions will contain read group IDs, and will end in `.bam`."
        }
    }

    parameter_meta {
        bam: {
            description: "Input BAM format file to split",
            stream: true,
        }
        prefix: "Prefix for the split BAM files. The extensions will contain read group IDs, and will end in `.bam`."
        reject_unaccounted_reads: {
            description: "If true, error if there are reads present that do not have read group information matching the header.",
            group: "Common",
        }
        reject_empty_output: {
            description: "If true, error if any output BAMs are empty.",
            group: "Common",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean reject_unaccounted_reads = true
        Boolean reject_empty_output = true
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools split \
            --threads "$n_cores" \
            -u "~{prefix}.unaccounted_reads.bam" \
            -f '~{prefix}.%!.bam' \
            "~{bam}"

        samtools head \
            --threads "$n_cores" \
            -h 0 \
            -n 1 \
            "~{prefix}.unaccounted_reads.bam" \
            > first_unaccounted_read.sam

        EXITCODE=0
        if ~{reject_unaccounted_reads} && [ -s first_unaccounted_read.sam ]; then
            >&2 echo "There are reads present with bad or missing RG tags!"
            EXITCODE=21
        else
            rm "~{prefix}.unaccounted_reads.bam"
        fi
        rm first_unaccounted_read.sam

        # Check that all output BAMs have at least
        # one read in them.
        if ~{reject_empty_output}; then
            for out_bam in *.bam; do
                samtools head \
                    --threads "$n_cores" \
                    -h 0 \
                    -n 1 \
                    "$out_bam" \
                    > first_read.sam

                if [ ! -s first_read.sam ]; then
                    >&2 echo "No reads are in output BAM $out_bam!"
                    >&2 echo "This is likely caused by malformed RG records."
                    EXITCODE=42
                fi
            done
            rm first_read.sam
        fi
        
        exit $EXITCODE
    >>>

    output {
       Array[File] split_bams = glob("*.bam")
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
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
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".flagstat.txt"
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools flagstat --threads "$n_cores" "~{bam}" > "~{outfile_name}"
    >>>

    output {
       File flagstat_report = outfile_name
    }

    runtime {
        memory: "5 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
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
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools index --threads "$n_cores" "~{bam}" "~{outfile_name}"
    >>>

    output {
        File bam_index = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task subsample {
    meta {
        description: "Randomly subsamples the input BAM, in order to produce an output BAM with approximately the desired number of reads."
        help: "A `desired_reads` **greater than zero** must be supplied. A `desired_reads <= 0` will result in task failure. Sampling is probabalistic and will be approximate to `desired_reads`. Read count will not be exact. A `sampled_bam` will not be produced if the input BAM read count is less than or equal to `desired_reads`."
        outputs: {
            orig_read_count: "A TSV report containing the original read count before subsampling. If subsampling was requested but the input BAM had less than `desired_reads`, no read count will be filled in (instead there will be a `dash`).",
            sampled_bam: "The subsampled input BAM. Only present if subsampling was performed.",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to subsample"
        desired_reads: "How many reads should be in the ouput BAM? Output BAM read count will be approximate to this value. **Must be greater than zero.** A `desired_reads <= 0` will result in task failure."
        prefix: "Prefix for the BAM file. The extension `.sampled.bam` will be added."
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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

    String suffixed = prefix + ".sampled"

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        if [[ ~{desired_reads} -le 0 ]]; then
            >&2 echo "'desired_reads' must be greater than zero!"
            >&2 echo "Task failed!"
            exit 42
        fi

        read_count="$(samtools view \
            --threads "$n_cores" \
            --count \
            "~{bam}"
        )"

        if [[ "$read_count" -eq 0 ]]; then
            >&2 echo "Read count is 0. Cannot subsample."
            >&2 echo "Please check that the input BAM is not empty."
            >&2 echo "Command failed!"
            exit 42
        fi

        if [[ "$read_count" -gt "~{desired_reads}" ]]; then
            # the BAM has more than ~{desired_reads} reads, meaning we should
            # subsample it.
            frac=$( \
                awk -v desired_reads=~{desired_reads} \
                    -v read_count="$read_count" \
                        'BEGIN{
                            printf "%1.8f",
                            ( desired_reads / read_count )
                        }' \
            )
            samtools view \
                --threads "$n_cores" \
                -hb \
                -s "$frac" \
                "~{bam}" \
                > "~{suffixed}.bam"

            # Report the original read count.
            # Use 'prefix' as the entry name.
            # Use the '.sampled' suffixed name in the filename
            # because that is the name of the output BAM.
            # This might seem odd but it works best with MultiQC.
            {
                echo -e "sample\toriginal read count"
                echo -e "~{prefix}\t$read_count"
            } > "~{suffixed}.orig_read_count.tsv"
        else
            # the BAM has less than or equal to ~{desired_reads} reads,
            # meaning we should just use it directly without subsampling.

            # Do not report an original read count,
            # as it is the same as the input BAM. Just write a dash.
            # Do not use the '.sampled' suffixed name in the filename.
            # Use ~{prefix} instead, so it matches the any downstream calls
            # which also use ~{prefix}.
            # This is for MultiQC purposes.
            {
                echo -e "sample\toriginal read count"
                echo -e "~{prefix}\t-"
            } > "~{prefix}.orig_read_count.tsv"
        fi

        # Check that if output was created,
        # it contains at least one read.
        if [ -e "~{suffixed}.bam" ]; then
            samtools head \
                --threads "$n_cores" \
                -h 0 \
                -n 1 \
                "~{suffixed}.bam" \
                > first_read.sam

            if [ ! -s first_read.sam ]; then
                >&2 echo "No reads are in the output BAM!"
                >&2 echo "This should not be possible! Please report this as a bug."
                rm first_read.sam
                exit 42
            fi
            rm first_read.sam
        fi

    >>>

    output {
        File orig_read_count = glob("*.orig_read_count.tsv")[0]
        File? sampled_bam = suffixed + ".bam"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task filter {
    meta {
        description: "Filters a BAM based on its bitwise flag value."
        help: "This task is a wrapper around `samtools view`. This task will fail if there are no reads in the output BAM. This can happen either because the input BAM was empty or because the supplied `bitwise_filter` was too strict. If you want to down-sample a BAM, use the `subsample` task instead."
        outputs: {
            filtered_bam: "BAM file that has been filtered based on the input flags"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to filter"
        bitwise_filter: "A set of 4 possible read filters to apply. This is a `FlagFilter` object (see ../data_structures/flag_filter.wdl for more information)."
        prefix: "Prefix for the filtered BAM file. The extension `.bam` will be added."
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        FlagFilter bitwise_filter
        String prefix = basename(bam, ".bam") + ".filtered"
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools view \
            --threads "$n_cores" \
            -hb \
            -f "~{bitwise_filter.include_if_all}" \
            -F "~{bitwise_filter.exclude_if_any}" \
            --rf "~{bitwise_filter.include_if_any}" \
            -G "~{bitwise_filter.exclude_if_all}" \
            "~{bam}" \
            > "~{prefix}.bam"

        samtools head \
            --threads "$n_cores" \
            -h 0 \
            -n 1 \
            "~{prefix}.bam" \
            > first_read.sam

        if [ ! -s first_read.sam ]; then
            >&2 echo "No reads are in the output BAM!"
            >&2 echo "Please check that the input BAM is not empty"
            >&2 echo "and that the filter is not too restrictive."
            >&2 echo "Command failed!"
            rm first_read.sam
            exit 42
        fi
        rm first_read.sam
    >>>

    output {
        File filtered_bam = prefix + ".bam"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
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
            group: "Common",
        }
        name_sorted: {
            description: "Are _all_ input BAMs `queryname` sorted (true)? Or are _all_ input BAMs `coordinate` sorted (false)?",
            group: "Common",
        }
        combine_rg: {
            description: "When several input files contain @RG headers with the same ID, emit only one of them (namely, the header line from the first file we find that ID in) to the merged output file. Combining these similar headers is usually the right thing to do when the files being merged originated from the same file. Without `-c`, all @RG headers appear in the output file, with random suffixes added to their IDs where necessary to differentiate them.",
            group: "Common",
        }
        combine_pg: {
            description: "Similarly to `combine_rg`: for each @PG ID in the set of files to merge, use the @PG line of the first file we find that ID in rather than adding a suffix to differentiate similar IDs.",
            group: "Common",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        bams=""
        for file in ~{sep(" ", squote(bams))}
        do
          # This will fail (intentionally) if there are duplicate names
          # in the input BAM array.
          ln -s "$file" .
          bams+=" $(basename "$file")"
        done

        samtools merge \
            --threads "$n_cores" \
            ~{"-h \"" + new_header + "\""} \
            ~{if name_sorted then "-n" else ""} \
            ~{if (region != "") then "-R \"" + region + "\"" else ""} \
            ~{if attach_rg then "-r" else ""} \
            ~{if combine_rg then "-c" else ""} \
            ~{if combine_pg then "-p" else ""} \
            "~{prefix}.bam" \
            "${bams[@]}"
    >>>

    output {
        File merged_bam = prefix + ".bam"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
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
            group: "Common",
        }
        prefix: "Prefix for the BAM file. The extension `.bam` will be added."
        orphan_only: {
            description: "Only add RG tags to orphans (true)? Or _also_ overwrite all existing RG tags (including any in the header) (false)?",
            group: "Common",
        }
        overwrite_header_record: {
            description: "Overwrite an existing @RG line, if a new one with the same ID value is provided?",
            group: "Common",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools addreplacerg \
            --threads "$n_cores" \
            ~{sep(" ", prefix("-r ", squote(read_group_line)))} \
            ~{"-R \"" + read_group_id + "\""} \
            -m ~{if orphan_only then "orphan_only" else "overwrite_all"} \
            ~{if overwrite_header_record then "-w" else ""} \
            -o "~{outfile_name}" \
            "~{bam}"
    >>>

    output {
        File tagged_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
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
            group: "Common",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools collate \
            --threads "$n_cores" \
            ~{if fast_mode then "-f" else ""} \
            -o "~{outfile_name}" \
            "~{bam}"
    >>>

    output {
        File collated_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task bam_to_fastq {
    meta {
        description: "Converts an input BAM file into FASTQ(s) using `samtools fastq`."
        help: "If `paired_end == false`, then _all_ reads in the BAM will be output to a single FASTQ file. Use `bitwise_filter` argument to remove any unwanted reads. An exit-code of `42` indicates that no reads were present in the output FASTQs. An exit-code of `43` indicates that unexpected reads were discovered in the input BAM."
        outputs: {
            collated_bam: "A collated BAM (reads sharing a name next to each other, no other guarantee of sort order). Only generated if `retain_collated_bam` and `paired_end` are both true. Has the name `~{prefix}.collated.bam`.",
            read_one_fastq_gz: "Gzipped FASTQ file with 1st reads in pair. Only generated if `paired_end` is true. Has the name `~{prefix}.R1.fastq.gz`.",
            read_two_fastq_gz: "Gzipped FASTQ file with 2nd reads in pair. Only generated if `paired_end` is true. Has the name `~{prefix}.R2.fastq.gz`.",
            singleton_reads_fastq_gz: "Gzipped FASTQ containing singleton reads. Only generated if `paired_end` and `output_singletons` are both true. Has the name `~{prefix}.singleton.fastq.gz`.",
            single_end_reads_fastq_gz: "A gzipped FASTQ containing all reads. Only generated if `paired_end` is false. Has the name `~{prefix}.fastq.gz`.",
        }
    }

    parameter_meta {
        bam: "Input BAM format file to convert to FASTQ(s)"
        bitwise_filter: "A set of 4 possible read filters to apply during conversion to FASTQ. This is a `FlagFilter` object (see ../data_structures/flag_filter.wdl for more information). By default, it will **remove secondary and supplementary reads** from the output FASTQs."
        prefix: "Prefix for the collated BAM and FASTQ files. The extensions `.collated.bam` and `[,.R1,.R2,.singleton].fastq.gz` will be added."
        paired_end: {
            description: "Is the data Paired-End? If `paired_end == false`, then _all_ reads in the BAM will be output to a single FASTQ file. Use `bitwise_filter` argument to remove any unwanted reads.",
            group: "Common",
        }
        collated: {
            description: "Is the BAM collated (or name-sorted)? If `collated == true`, then the input BAM will be run through `samtools fastq` without preprocessing. If `collated == false`, then `samtools collate` must be run on the input BAM before conversion to FASTQ. Ignored if `paired_end == false`.",
            group: "Common",
        }
        retain_collated_bam: {
            description: "Save the collated BAM to disk and output it (true)? This slows performance and **substantially** increases storage requirements. Be aware that collated BAMs occupy much more space than either position sorted or name sorted BAMs (due to the compression algorithm). Ignored if `collated == true` **or** `paired_end == false`.",
            group: "Common",
        }
        fast_mode: {
            description: "Fast mode for `samtools collate`? If `true`, this **removes secondary and supplementary reads** during the `collate` step. If `false`, secondary and supplementary reads will be retained in the `collated_bam` output (if created). Defaults to the opposite of `retain_collated_bam`. Ignored if `collated == true` **or** `paired_end == false`.",
            group: "Common",
        }
        append_read_number: {
            description: "Append /1 and /2 suffixes to read names?",
            group: "Common",
        }
        output_singletons: "Output singleton reads as their own FASTQ? Ignored if `paired_end == false`."
        fail_on_unexpected_reads: {
            description: "Should the task fail if reads with an unexpected `first`/`last` bit setting are discovered?",
            help: "The definition of 'unexpected' depends on whether the values of `paired_end` and `output_singletons` are true or false. If `paired_end` is `false`, no reads are considered unexpected, and _every_ read (not caught by `bitwise_filter`) will be present in the resulting FASTQ regardless of `first`/`last` bit settings. This setting will be ignored in that case. If `paired_end` is `true` then reads that don't satisfy `first` XOR `last` are considered unexpected (i.e. reads that have neither `first` nor `last` set or reads that have both `first` and `last` set). If `output_singletons` is `false`, singleton reads are considered unexpected. A singleton read is a read with either the `first` or the `last` bit set (but not both) and that possesses a _unique_ QNAME; i.e. it is a read without a pair when all reads are expected to be paired. But if `output_singletons` is `true`, these singleton reads will be output as their own FASTQ instead of causing the task to fail. If `fail_on_unexpected_reads` is `false`, then all the above cases will be ignored. Any 'unexpected' reads will be silently discarded.",
            group: "Common",
        }
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        FlagFilter bitwise_filter = FlagFilter {
            include_if_all: "0x0",
            exclude_if_any: "0x900",
            include_if_any: "0x0",
            exclude_if_all: "0x0",
        }
        String prefix = basename(bam, ".bam")
        Boolean paired_end = true
        Boolean collated = false
        Boolean retain_collated_bam = false
        Boolean fast_mode = !retain_collated_bam
        Boolean append_read_number = true
        Boolean output_singletons = false
        Boolean fail_on_unexpected_reads = false
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = (
        if (collated || !paired_end)
        then 4
        else (ceil(bam_size * 0.4) + 4)
    ) + modify_memory_gb
    Int disk_size_gb = ceil(bam_size * (
        if (retain_collated_bam && !collated && paired_end)
        then 5
        else 2
    )) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        mkfifo bam_pipe
        if ! ~{collated} && ~{paired_end}; then
            samtools collate \
                ~{if retain_collated_bam then "" else "-u"} \
                --threads "$n_cores" \
                ~{if fast_mode then "-f" else ""} \
                -O \
                "~{bam}" \
                | tee ~{(
                    if retain_collated_bam
                    then "\"" + prefix + ".collated.bam\""
                    else ""
                )} \
                > bam_pipe \
                &
        else
            samtools view -h --threads "$n_cores" "~{bam}" > bam_pipe &
        fi

        samtools fastq \
            --threads "$n_cores" \
            -f "~{bitwise_filter.include_if_all}" \
            -F "~{bitwise_filter.exclude_if_any}" \
            --rf "~{bitwise_filter.include_if_any}" \
            -G "~{bitwise_filter.exclude_if_all}" \
            ~{(
                if append_read_number
                then "-N"
                else "-n"
            )} \
            -1 ~{(
                if paired_end
                then "\"" + prefix + ".R1.fastq.gz\""
                else "\"" + prefix + ".fastq.gz\""
            )} \
            -2 ~{(
                if paired_end
                then "\"" + prefix + ".R2.fastq.gz\""
                else "\"" + prefix + ".fastq.gz\""
            )} \
            ~{(
                if paired_end
                then (
                    if output_singletons
                    then "-s \"" + prefix + ".singleton.fastq.gz\""
                    else "-s junk.singleton.fastq.gz"
                )
                else ""
            )} \
            -0 ~{(
                if paired_end
                then "junk.unknown_bit_setting.fastq.gz"
                else "\"" + prefix + ".fastq.gz\""
            )} \
            bam_pipe

        rm bam_pipe

        # Check that some output is non-empty
        if [ -z "$(gunzip -c "~{prefix}"*.fastq.gz | head -c 1 | tr '\0\n' __)" ]; then
            >&2 echo "No reads are in any output FASTQ"
            >&2 echo "Command failed!"
            exit 42
        fi

        # Check that there weren't any unexpected reads in the input BAM
        if ~{fail_on_unexpected_reads} \
            && [ -n "$(gunzip -c junk.*.fastq.gz | head -c 1 | tr '\0\n' __)" ]
        then
            >&2 echo "Discovered unexpected reads"
            exit 43
        fi
    >>>

    output {
        File? collated_bam = "~{prefix}.collated.bam"
        File? read_one_fastq_gz = "~{prefix}.R1.fastq.gz"
        File? read_two_fastq_gz = "~{prefix}.R2.fastq.gz"
        File? singleton_reads_fastq_gz = "~{prefix}.singleton.fastq.gz"
        File? single_end_reads_fastq_gz = "~{prefix}.fastq.gz"
    }

    runtime {
        cpu: ncpu
        disks: "~{disk_size_gb} GB"
        memory: "~{memory_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task fixmate {
    meta {
        description: "Runs `samtools fixmate` on the name-collated input BAM file. This fills in mate coordinates and insert size fields among other tags and fields."
        help: "This task assumes a name-sorted or name-collated input BAM. If you have a position-sorted BAM, please use the `position_sorted_fixmate` task. This task runs `fixmate` and outputs a BAM in the same order as the input."
        outputs: {
            fixmate_bam: "The BAM resulting from running `samtools fixmate` on the input BAM"
        }
    }

    parameter_meta {
        bam: {
            description: "Input BAM format file to add mate information. Must be name-sorted or name-collated.",
            stream: true,
        }
        prefix: "Prefix for the output file. The extension specified with the `extension` parameter will be added."
        extension: {
            description: "File format extension to use for output file.",
            choices: [
                ".bam",
                ".cram",
            ],
            group: "Common",
        }
        add_cigar: {
            description: "Add template cigar `ct` tag",
            tool_default: false,
            group: "Common",
        }
        add_mate_score: {
            description: "Add mate score tags. These are used by `markdup` to select the best reads to keep.",
            tool_default: false,
            group: "Common",
        }
        disable_flag_sanitization: "Disable all flag sanitization?"
        disable_proper_pair_check: "Disable proper pair check [ensure one forward and one reverse read in each pair]"
        remove_unaligned_and_secondary: "Remove unmapped and secondary reads"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".fixmate"
        String extension = ".bam"
        Boolean add_cigar = true
        Boolean add_mate_score = true
        Boolean disable_flag_sanitization = false
        Boolean disable_proper_pair_check = false
        Boolean remove_unaligned_and_secondary = false
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
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools fixmate \
            --threads "$n_cores" \
            ~{if remove_unaligned_and_secondary then "-r" else ""} \
            ~{if disable_proper_pair_check then "-p" else ""} \
            ~{if add_cigar then "-c" else ""} \
            ~{if add_mate_score then "-m" else ""} \
            ~{if disable_flag_sanitization then "-z off" else ""} \
            "~{bam}" \
            "~{prefix}~{extension}"
    >>>

    output {
        File fixmate_bam = "~{prefix}~{extension}"
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task position_sorted_fixmate {
    meta {
        description: "Runs `samtools fixmate` on the position-sorted input BAM file and output a position-sorted BAM. `fixmate` fills in mate coordinates and insert size fields among other tags and fields. `samtools fixmate` assumes a name-sorted or name-collated input BAM. If you already have a collated BAM, please use the `fixmate` task. This task collates the input BAM, runs `fixmate`, and then resorts the output into a position-sorted BAM."
        outputs: {
            fixmate_bam: "BAM file with mate information added"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to add mate information. Must be position-sorted."
        prefix: "Prefix for the output file. The extension `.bam` will be added."
        fast_mode: {
            description: "Use fast mode (output primary alignments only)?",
            group: "Common",
        }
        add_cigar: {
            description: "Add template cigar `ct` tag",
            tool_default: false,
            group: "Common",
        }
        add_mate_score: {
            description: "Add mate score tags. These are used by `markdup` to select the best reads to keep.",
            tool_default: false,
            group: "Common",
        }
        disable_flag_sanitization: "Disable all flag sanitization?"
        disable_proper_pair_check: "Disable proper pair check [ensure one forward and one reverse read in each pair]?"
        remove_unaligned_and_secondary: "Remove unmapped and secondary reads"
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".fixmate"
        Boolean fast_mode = false
        Boolean add_cigar = true
        Boolean add_mate_score = true
        Boolean disable_flag_sanitization = false
        Boolean disable_proper_pair_check = false
        Boolean remove_unaligned_and_secondary = false
        Boolean use_all_cores = false
        Int ncpu = 2
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = ceil(bam_size * 0.2) + 4 + modify_memory_gb
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools collate \
            --threads "$n_cores" \
            ~{if fast_mode then "-f" else ""} \
            -u \
            -O \
            "~{bam}" \
            | samtools fixmate \
                --threads "$n_cores" \
                -u \
                ~{if remove_unaligned_and_secondary then "-r" else ""} \
                ~{if disable_proper_pair_check then "-p" else ""} \
                ~{if add_cigar then "-c" else ""} \
                ~{if add_mate_score then "-m" else ""} \
                ~{if disable_flag_sanitization then "-z off" else ""} \
                - \
                - \
                | samtools sort \
                    --threads "$n_cores" \
                    -o "~{prefix + ".bam"}" \
                    -
    >>>

    output {
        File fixmate_bam = "~{prefix}.bam"
    }

    runtime {
        cpu: ncpu
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

#@ except: MatchingOutputMeta
task markdup {
    meta {
        description: "**[DEPRECATED]** Runs `samtools markdup` on the position-sorted input BAM file. This creates a report and optionally a new BAM with duplicate reads marked."
        help: "This task assumes `samtools fixmate` has already been run on the input BAM. If it has not, then the output may be incorrect. A name-sorted or collated BAM can be run through the `fixmate` task (and then position-sorted prior to this task) or a position-sorted BAM can be run through the `position_sorted_fixmate` task. Deprecated due to extremely high memory usage for certain RNA-Seq samples when searching for optical duplicates. Use `mark_duplicates` in `./picard.wdl` instead."
        deprecated: true
    }

    parameter_meta {
        bam: "Input BAM format file to mark duplicates in"
        prefix: "Prefix for the output file."
        read_coords_regex: {
            description: "Regular expression to extract read coordinates from the QNAME field. This takes a POSIX regular expression for at least x and y to be used in optical duplicate marking It can also include another part of the read name to test for equality, eg lane:tile elements. Elements wanted are captured with parentheses. The default is meant to capture information from Illumina style read names. Ignored if `optical_distance == 0`. If changing `read_coords_regex`, make sure that `coordinates_order` matches.",
            tool_default: "`([!-9;-?A-~]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+):([0-9]+):([0-9]+)`",
        }
        coordinates_order: {
            description: "The order of the elements captured in the `read_coords_regex` regular expression. Default is `txy` where `t` is a part of the read name selected for string comparison and `x`/`y` are the coordinates used for optical duplicate detection. Ignored if `optical_distance == 0`.",
            choices: [
                "txy",
                "tyx",
                "xyt",
                "yxt",
                "xty",
                "ytx",
                "xy",
                "yx",
            ],
        }
        create_bam: "Create a new BAM with duplicate reads marked? If `false`, then only a markdup report will be generated."
        remove_duplicates: "Remove duplicates from the output BAM? Ignored if `create_bam == false`."
        mark_supp_or_sec_or_unmapped_as_duplicates: "Mark supplementary, secondary, or unmapped alignments of duplicates as duplicates? As this takes a quick second pass over the data it will increase running time. Ignored if `create_bam == false`."
        json: "Output a JSON report instead of a text report? Either are parseable by MultiQC."
        mark_duplicates_with_do_tag: "Mark duplicates with the `do` (`d`uplicate `o`riginal) tag? The `do` tag contains the name of the \"original\" read that was duplicated. Ignored if `create_bam == false`."
        duplicate_count: "Record the original primary read duplication count (include itself) in a `dc` tag? Ignored if `create_bam == false`."
        include_qc_fails: "Include reads that have the QC-failed flag set in duplicate marking? This can increase the number of duplicates found. Ignored if `create_bam == false`."
        duplicates_of_duplicates_check: "Check duplicates of duplicates for correctness? Performs further checks to make sure all optical duplicates are found. Also operates on `mark_duplicates_with_do_tag` tagging where reads may be tagged with the best quality read. Disabling this option can speed up duplicate marking when there are a great many duplicates for each original read. Ignored if `create_bam == false` **or** `optical_distance == 0`."
        use_read_groups: "Only mark duplicates _within_ the same Read Group? Ignored if `create_bam == false`."
        use_all_cores: {
            description: "Use all cores? Recommended for cloud environments.",
            group: "Common",
        }
        max_readlen: "Expected maximum read length."
        optical_distance: "Maximum distance between read coordinates to consider them optical duplicates. If `0`, then optical duplicate marking is disabled. Suggested settings of 100 for HiSeq style platforms or about 2500 for NovaSeq ones. When set above `0`, duplicate reads are tagged with `dt:Z:SQ` for optical duplicates and `dt:Z:LB` otherwise. Calculation of distance depends on coordinate data embedded in the read names, typically produced by the Illumina sequencing machines. Optical duplicate detection will not work on non-standard names without modifying `read_coords_regex`. If changing `read_coords_regex`, make sure that `coordinates_order` matches."
        ncpu: {
            description: "Number of cores to allocate for task",
            group: "Common",
        }
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam") + ".markdup"
        String read_coords_regex = "[!-9;-?A-~:]+:([!-9;-?A-~]+):([0-9]+):([0-9]+)"
        String coordinates_order = "txy"
        Boolean create_bam = true
        Boolean remove_duplicates = false
        Boolean mark_supp_or_sec_or_unmapped_as_duplicates = false
        Boolean json = false
        Boolean mark_duplicates_with_do_tag = false
        Boolean duplicate_count = false
        Boolean include_qc_fails = false
        Boolean duplicates_of_duplicates_check = false  # tool default is true
        Boolean use_read_groups = false
        Boolean use_all_cores = false
        Int max_readlen = 300
        Int optical_distance = 0
        Int ncpu = 2
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = ceil(bam_size * 3) + 4 + modify_memory_gb
    Int disk_size_gb = ceil(bam_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi
        # -1 because samtools uses one more core than `--threads` specifies
        (( n_cores -= 1 ))

        samtools markdup \
            --threads "$n_cores" \
            -f "~{prefix + if json then ".json" else ".txt"}" \
            --read-coords '~{read_coords_regex}' \
            --coords-order "~{coordinates_order}" \
            ~{if remove_duplicates then "-r" else ""} \
            ~{if mark_supp_or_sec_or_unmapped_as_duplicates then "-S" else ""} \
            ~{if mark_duplicates_with_do_tag then "-t" else ""} \
            ~{if duplicate_count then "--duplicate-count" else ""} \
            ~{if include_qc_fails then "--include-fails" else ""} \
            ~{if duplicates_of_duplicates_check then "" else "--no-multi-dup"} \
            ~{if use_read_groups then "--use-read-groups" else ""} \
            -l ~{max_readlen} \
            -d ~{optical_distance} \
            -c \
            "~{bam}" \
            "~{if create_bam then prefix + ".bam" else "/dev/null"}"
    >>>

    output {
        File markdup_report = prefix + if json then ".json" else ".txt"
        File? markdup_bam = prefix + ".bam"
    }

    runtime {
        cpu: ncpu
        disks: "~{disk_size_gb} GB"
        memory: "~{memory_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task faidx {
    meta {
        description: "Creates a `.fai` FASTA index for the input FASTA"
        outputs: {
            fasta_index: "A `.fai` FASTA index associated with the input FASTA. Filename will be `basename(fasta) + '.fai'`."
        }
    }

    parameter_meta {
        fasta: "Input FASTA format file to index. Optionally gzip compressed."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File fasta
        Int modify_disk_size_gb = 0
    }

    Float fasta_size = size(fasta, "GiB")
    Int disk_size_gb = ceil(fasta_size * 2.5) + 10 + modify_disk_size_gb

    String outfile_name = basename(fasta, ".gz") + ".fai"

    command <<<
        set -euo pipefail

        ref_fasta=~{basename(fasta, ".gz")}
        gunzip -c "~{fasta}" > "$ref_fasta" \
            || ln -sf "~{fasta}" "$ref_fasta"

        samtools faidx -o "~{outfile_name}" "$ref_fasta"
    >>>

    output {
        File fasta_index = outfile_name
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}
