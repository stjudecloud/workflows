## # SAMtools
##
## This WDL tool wraps the [SAMtools package](http://samtools.sourceforge.net/).
## SAMtools provides utlities for manipulating SAM format sequence alignments.

version 1.0

task quickcheck {
    input {
        File bam
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        samtools quickcheck ~{bam}
    }

    runtime {
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
        maxRetries: max_retries
    }

    output {
        File checked_bam = bam
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool runs Samtools quickcheck on the input BAM file. This checks that the BAM file appears to be intact, e.g. header exists, at least one sequence is present, and the end-of-file marker exists."
    }

    parameter_meta {
        bam: "Input BAM format file to quickcheck"
    }
}

task split {
    input {
        File bam
        Int ncpu = 1
        Boolean reject_unaccounted = true
        String prefix = basename(bam, ".bam")
        Int max_retries = 1
        Int? disk_size_gb
        Boolean use_all_cores = false
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])

    command {
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        samtools split --threads "$n_cores" -u ~{prefix}.unaccounted_reads.bam -f '%*_%!.%.' ~{bam}
        samtools view  --threads "$n_cores" ~{prefix}.unaccounted_reads.bam > unaccounted_reads.bam
        if ~{reject_unaccounted} && [ -s unaccounted_reads.bam ]
            then exit 1; 
            else rm ~{prefix}.unaccounted_reads.bam
        fi 
        rm unaccounted_reads.bam
    }
 
    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
        maxRetries: max_retries
    }

    output {
       Array[File] split_bams = glob("*.bam")
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool runs Samtools split on the input BAM file. This splits the BAM by read group into one or more output files. It optionally errors if there are reads present that do not belong to a read group."
    }

    parameter_meta {
        bam: "Input BAM format file to split"
        reject_unaccounted: "If true, error if there are reads present that do not have read group information."
    }
}

task flagstat {
    input {
        File bam
        String outfile_name = basename(bam, ".bam")+".flagstat.txt"
        Int max_retries = 1
        Int memory_gb = 5
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        samtools flagstat ~{bam} > ~{outfile_name}
    }

    runtime {
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
        memory: memory_gb + " GB"
        maxRetries: max_retries
    }

    output { 
       File flagstat_report = outfile_name
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates a `samtools flagstat` report for the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate flagstat for"
    }
}

task index {
    input {
        File bam
        Int ncpu = 1
        String outfile_name = basename(bam)+".bai"
        Int max_retries = 1
        Int memory_gb = 15
        Boolean use_all_cores = false
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command {
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        samtools index -@ "$n_cores" ~{bam} ~{outfile_name}
    }

    runtime {
        cpu: ncpu
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
        memory: memory_gb + " GB"
        dx_instance_type: "azure:mem2_ssd1_x4"
        maxRetries: max_retries
    }

    output {
        File bam_index = outfile_name
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool runs Samtools flagstat on the input BAM file. Produces statistics about the alignments based on the bit flags set in the BAM."
    }

    parameter_meta {
        bam: "Input BAM format file to index"
    }
}

task subsample {
    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".subsampled.bam"
        Int desired_reads = 500000
        Int max_retries = 1
        Int ncpu = 1
        Boolean use_all_cores = false
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        if [[ "$(samtools view --threads "$n_cores" ~{bam} | head -n ~{desired_reads} | wc -l)" -ge "~{desired_reads}" ]]; then
            # the BAM has at least ~{desired_reads} reads, meaning we should
            # subsample it.
            initial_frac=0.00001
            initial_reads=$(samtools view --threads "$n_cores" -s "$initial_frac" ~{bam} | wc -l)
            frac=$( \
                awk -v desired_reads=~{desired_reads} \
                    -v initial_reads="$initial_reads" \
                    -v initial_frac="$initial_frac" \
                        'BEGIN{printf "%1.8f", ( desired_reads / initial_reads * initial_frac )}' \
                )
            samtools view --threads "$n_cores" -h -b -s "$frac" ~{bam} > ~{outfile_name}
        else
            # the BAM has less than ~{desired_reads} reads, meaning we should
            # just use it directly without subsampling.
            mv ~{bam} ~{outfile_name}
        fi
    >>>

    output {
        File sampled_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
        maxRetries: max_retries
    }
}

task merge {
    input {
        Array[File] bams
        String outfile_name = basename(bams[0], ".bam") + ".merged.bam"
        File? new_header
        Boolean attach_rg = true
        Int max_retries = 1
        Int ncpu = 1
        Boolean use_all_cores = false
    }

    Float bam_size = size(bams, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)
    String rg_arg = if attach_rg then "-r" else ""

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        header_arg=""
        if [ -e ~{new_header} ]
        then
            header_arg="-h ~{new_header}"
        fi

        samtools merge --threads "$n_cores" $header_arg ~{rg_arg} ~{outfile_name} ~{sep=' ' bams}

    >>>

    output {
        File merged_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
        maxRetries: max_retries
    }
}

task addreplacerg {
    input {
        File bam
        String outfile_name = basename(bam, ".bam") + ".read_group.bam"
        String read_group_id
        Int max_retries = 1
        Int ncpu = 1
        Boolean use_all_cores = false
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        samtools addreplacerg --threads "$n_cores" -R ~{read_group_id} -o ~{outfile_name} ~{bam}
    >>>

    output {
        File tagged_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disk: disk_size + " GB"
        docker: 'ghcr.io/stjudecloud/samtools:1.0.2'
        maxRetries: max_retries
    }


    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
        description: "This WDL tool runs Samtools addreplacerg on the input BAM file. This adds an existing read group record to reads in the BAM lacking read group tags."
    }

    parameter_meta {
        bam: "Input BAM format file to add read group information"
        read_group_id: "Existing read group ID in BAM to add to reads"
    }
}

task collate {
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean use_all_cores = false
        Int memory_gb = 24
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    parameter_meta {
        bam: "Input BAM format file to collate"
    }

    String outfile_name = prefix + ".collated.bam"

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        samtools collate --threads "$n_cores" -o ~{outfile_name} ~{bam}
    >>>

    output {
        File collated_bam = outfile_name
    }

    runtime {
        cpu: ncpu
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: max_retries
    }


    meta {
        description: "This WDL tool runs `samtools collate` on the input BAM file. Shuffles and groups reads together by their names."
    }
}

task bam_to_fastq {
    # TODO: enable SE mode and expose more parameters
    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean use_all_cores = false
        Int memory_gb = 4
        Int modify_disk_size_gb = 0
        Int ncpu = 1
        Int max_retries = 1
    }

    parameter_meta {
        bam: "Input BAM format file to split into read1 and read2 FastQs"
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil((bam_size * 2) + 10) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{use_all_cores}" = "true" ]; then
            n_cores=$(nproc)
        fi

        samtools fastq \
            --threads "$n_cores" \
            -1 ~{prefix}_R1.fastq.gz \
            -2 ~{prefix}_R2.fastq.gz \
            ~{bam}
    >>>

    output {
        File read1 = "~{prefix}_R1.fastq.gz"
        File read2 = "~{prefix}_R2.fastq.gz"
    }

    runtime {
        cpu: ncpu
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
        maxRetries: max_retries
    }


    meta {
        description: "This WDL tool runs `samtools fastq` on the input BAM file. Splits the BAM into read1 and read2 FastQ files."
    }
}
