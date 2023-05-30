## # SAMtools
##
## This WDL file wraps the [SAMtools package](http://samtools.sourceforge.net/).
## SAMtools provides utlities for manipulating SAM format sequence alignments.

version 1.0

task quickcheck {
    meta {
        description: "This WDL task runs Samtools quickcheck on the input BAM file. This checks that the BAM file appears to be intact, e.g. header exists, at least one sequence is present, and the end-of-file marker exists."
    }

    parameter_meta {
        bam: "Input BAM format file to quickcheck"
    }

    input {
        File bam
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 1.2) + 10) + modify_disk_size_gb

    command <<<
        samtools quickcheck ~{bam}
    >>>

    output {
        File checked_bam = bam
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task split {
    meta {
        description: "This WDL task runs Samtools split on the input BAM file. This splits the BAM by read group into one or more output files. It optionally errors if there are reads present that do not belong to a read group."
    }

    parameter_meta {
        bam: "Input BAM format file to split"
        reject_unaccounted: "If true, error if there are reads present that do not have read group information."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean reject_unaccounted = true
        Boolean use_all_cores = false
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools split \
            --threads "$n_cores" \
            -u ~{prefix}.unaccounted_reads.bam \
            -f '%*_%!.%.' \
            ~{bam}
        
        samtools head \
            --threads "$n_cores" \
            -h 0 \
            -n 1 \
            ~{prefix}.unaccounted_reads.bam \
            > unaccounted_reads.bam
        
        if ~{reject_unaccounted} && [ -s unaccounted_reads.bam ]; then
            exit 1
        else
            rm ~{prefix}.unaccounted_reads.bam
        fi 
        rm unaccounted_reads.bam
    >>>

    output {
       Array[File] split_bams = glob("*.bam")
    }
 
    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task flagstat {
    meta {
        description: "This WDL tool generates a `samtools flagstat` report for the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate flagstat for"
    }

    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".flagstat.txt"
        Int memory_gb = 5
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 1.2) + 10) + modify_disk_size_gb

    command <<<
        samtools flagstat ~{bam} > ~{outfile_name}
    >>>

    output { 
       File flagstat_report = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task index {
    meta {
        description: "This WDL task runs Samtools flagstat on the input BAM file. Produces statistics about the alignments based on the bit flags set in the BAM."
    }

    parameter_meta {
        bam: "Input BAM format file to index"
    }

    input {
        File bam
        String outfile_name = basename(bam) + ".bai"  # TODO we might want to disable customizing the name for the BAI? Breaks our best practices, but also please do not change the BAI name
        Boolean use_all_cores = false
        Int memory_gb = 15  # TODO this seems like more RAM than should be necessary
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 1.5) + 10) + modify_disk_size_gb

    command {
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools index -@ "$n_cores" ~{bam} ~{outfile_name}
    }

    output {
        File bam_index = outfile_name
    }

    runtime {
        dx_instance_type: "azure:mem2_ssd1_x4"
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task subsample {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean use_all_cores = false
        Int desired_reads = 500000
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        if [[ \
            "$(samtools head --threads "$n_cores" -h 0 -n ~{desired_reads} ~{bam} \
                | wc -l)" \
            -ge "~{desired_reads}" \
        ]]; then
            # the BAM has at least ~{desired_reads} reads, meaning we should
            # subsample it.
            initial_frac=0.00001
            initial_reads=$( \
                samtools view --threads "$n_cores" -s "$initial_frac" ~{bam} | wc -l \
            )
            frac=$( \
                awk -v desired_reads=~{desired_reads} \
                    -v initial_reads="$initial_reads" \
                    -v initial_frac="$initial_frac" \
                        'BEGIN{
                            printf "%1.8f",
                            ( desired_reads / initial_reads * initial_frac )
                        }' \
                )
            samtools view --threads "$n_cores" -hb -s "$frac" ~{bam} \
                > ~{prefix}.subsampled.bam
        else
            # the BAM has less than ~{desired_reads} reads, meaning we should
            # just use it directly without subsampling.
            true
        fi
        touch success
    >>>

    output {
        File success = "success"
        File? sampled_bam = prefix + ".subsampled.bam"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task merge {
    input {
        Array[File] bams
        File? new_header
        String outfile_name = basename(bams[0], ".bam") + ".merged.bam"  # TODO is this desired behavior?
        Boolean attach_rg = true
        Boolean use_all_cores = false
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bams_size = size(bams, "GiB")
    Int disk_size_gb = ceil((bams_size * 2) + 10) + modify_disk_size_gb
    
    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools merge \
            --threads "$n_cores" \
            ~{if defined(new_header) then "-h " + new_header else ""} \
            ~{if attach_rg then "-r" else ""} \
            ~{outfile_name} \
            ~{sep=' ' bams}
    >>>

    output {
        File merged_bam = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task addreplacerg {
    meta {
        description: "This WDL task runs Samtools addreplacerg on the input BAM file. This adds an existing read group record to reads in the BAM lacking read group tags."
    }

    parameter_meta {
        bam: "Input BAM format file to add read group information"
        read_group_id: "Existing read group ID in BAM to add to reads"
    }

    input {
        File bam
        String read_group_id
        String outfile_name = basename(bam, ".bam") + ".read_group.bam"
        Boolean use_all_cores = false
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools addreplacerg \
            --threads "$n_cores" \
            -R ~{read_group_id} \
            -o ~{outfile_name} \
            ~{bam}
    >>>

    output {
        File tagged_bam = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'
        maxRetries: max_retries
    }
}

task collate {
    meta {
        description: "This WDL task runs `samtools collate` on the input BAM file. Shuffles and groups reads together by their names."
        outputs: {
            collated_bam: "A collated BAM (reads sharing a name next to each other, no other guarantee of sort order)"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to collate"
        prefix: "Prefix for the collated BAM file. The extension `.collated.bam` will be added."
        f: "Fast mode (primary alignments only)"
        use_all_cores: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".collated.bam"
        Boolean f = true
        Boolean use_all_cores = false
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = ceil(bam_size * 0.2) + 4 + modify_memory_gb
    Int disk_size_gb = ceil((bam_size * 4)) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if ~{use_all_cores}; then
            n_cores=$(nproc)
        fi

        samtools collate \
            --threads "$n_cores" \
            ~{if f then "-f" else ""} \
            -o ~{outfile_name} \
            ~{bam}
    >>>

    output {
        File collated_bam = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: max_retries
    }
}

task bam_to_fastq {
    meta {
        description: "This WDL task runs `samtools fastq` on the input BAM file. Splits the BAM into FastQ files. Assumes either a name sorted or collated BAM. For splitting a position sorted BAM see `collate_to_fastq`."
        outputs: {
            read_one_fastq_gz: "Gzipped FastQ file with 1st reads in pair"
            read_two_fastq_gz: "Gzipped FastQ file with 2nd reads in pair"
            singleton_reads_fastq_gz: "A gzipped FastQ containing singleton reads"
            interleaved_reads_fastq_gz: "An interleaved gzipped paired-end FastQ"
        }
    }

    parameter_meta {
        bam: "Input name sorted or collated BAM format file to convert into FastQ(s)"
        prefix: "Prefix for output FastQ(s). Extensions `[,_R1,_R2,.singleton].fastq.gz` will be added depending on other options."
        f: "Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        F: "Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/). This defaults to 0x900 representing filtering of secondary and supplementary alignments."
        G: "Only EXCLUDE reads with all of the bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        paired_end: "Is the data paired-end?"
        interleaved: "Create an interleaved FastQ file from paired-end data?"
        output_singletons: "Output singleton reads as their own FastQ?"
        memory_gb: "RAM to allocate for task, specified in GB"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
        use_all_cores: "Use all available cores? Recommended for cloud environments. Not recommended for cluster environments."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        String f = "0"
        String F = "0x900"
        String G = "0"
        Boolean paired_end = true
        Boolean interleaved = false
        Boolean output_singletons = false
        Boolean use_all_cores = false
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

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
            -1 ~{if interleaved
                then prefix + ".fastq.gz"
                else prefix + "_R1.fastq.gz"
            } \
            -2 ~{
                if paired_end then (
                    if interleaved then prefix + ".fastq.gz" else prefix + "_R2.fastq.gz"
                )
                else "/dev/null"
            } \
            -s ~{
                if output_singletons
                then prefix+".singleton.fastq.gz"
                else "/dev/null"
            } \
            -0 /dev/null \
            ~{bam}
    >>>

    output {
        # one of `read_one_fastq_gz` or `interleaved_reads_fastq_gz` is
        # guaranteed to exist at the end of execution
        File? read_one_fastq_gz = "~{prefix}_R1.fastq.gz"
        File? read_two_fastq_gz = "~{prefix}_R2.fastq.gz"
        File? singleton_reads_fastq_gz = "~{prefix}.singleton.fastq.gz"
        File? interleaved_reads_fastq_gz = "~{prefix}.fastq.gz"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: max_retries
    }
}

task collate_to_fastq {
    meta {
        description: "This WDL task runs `samtools collate` on the input BAM file then converts it into FastQ(s) using `samtools fastq`."
        outputs: {
            collated_bam: "A collated BAM (reads sharing a name next to each other, no other guarantee of sort order)"
            read_one_fastq_gz: "Gzipped FastQ file with 1st reads in pair"
	        read_two_fastq_gz: "Gzipped FastQ file with 2nd reads in pair"
            singleton_reads_fastq_gz: "Gzipped FastQ containing singleton reads"
            interleaved_reads_fastq_gz: "Interleaved gzipped paired-end FastQ"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to collate and convert to FastQ(s)"
        prefix: "Prefix for the collated BAM and FastQ files. The extensions `.collated.bam` and `[,_R1,_R2,.singleton].fastq.gz` will be added."
        f: "Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        F: "Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/). This defaults to 0x900 representing filtering of secondary and supplementary alignments."
        G: "Only EXCLUDE reads with all of the bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x` (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0` (i.e. /^0[0-7]+/)."
        fast_mode: "Fast mode for `samtools collate` (primary alignments only)"
        store_collated_bam: "Save the collated BAM (true) or delete it after FastQ split (false)?"
        paired_end: "Is the data paired-end?"
        interleaved: "Create an interleaved FastQ file from paired-end data?"
        output_singletons: "Output singleton reads as their own FastQ?"
        use_all_cores: "Use all cores? Recommended for cloud environments. Not recommended for cluster environments."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        String f = "0"
        String F = "0x900"
        String G = "0"
        Boolean fast_mode = true
        Boolean store_collated_bam = false
        Boolean paired_end = true
        Boolean interleaved = false
        Boolean output_singletons = false
        Boolean use_all_cores = false
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int memory_gb = ceil(bam_size * 0.2) + 4 + modify_memory_gb
    Int disk_size_gb = ceil((bam_size * 5) + 10) + modify_disk_size_gb

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
                -1 ~{if interleaved
                    then prefix + ".fastq.gz"
                    else prefix + "_R1.fastq.gz"
                } \
                -2 ~{
                    if paired_end then (
                        if interleaved
                        then prefix + ".fastq.gz"
                        else prefix + "_R2.fastq.gz"
                    )
                    else "/dev/null"
                } \
                -s ~{
                    if output_singletons
                    then prefix + ".singleton.fastq.gz"
                    else "/dev/null"
                } \
                -0 /dev/null
    >>>

    output {
        # one of `read_one_fastq_gz` or `interleaved_reads_fastq_gz` is
        # guaranteed to exist at the end of execution
        File? collated_bam = "~{prefix}.collated.bam"
        File? read_one_fastq_gz = "~{prefix}_R1.fastq.gz"
        File? read_two_fastq_gz = "~{prefix}_R2.fastq.gz"
        File? singleton_reads_fastq_gz = "~{prefix}.singleton.fastq.gz"
        File? interleaved_reads_fastq_gz = "~{prefix}.fastq.gz"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: max_retries
    }
}
