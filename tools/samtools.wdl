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

    command <<<
        samtools quickcheck ~{bam}
    >>>

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/samtools:1.0.0'
        maxRetries: max_retries
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
        Boolean? reject_unaccounted
        String prefix = basename(bam, ".bam")
        Int max_retries = 1
        Int? disk_size_gb
        Boolean detect_nproc = false
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])

    command <<<
        set -euo pipefail
        
        n_cores=~{ncpu}
        if [ ${true='true' false='' detect_nproc} ]
        then
            n_cores=`nproc`
        fi

        samtools split --threads $n_cores -u ~{prefix}.unaccounted_reads.bam -f '%*_%!.%.' ~{bam}
        samtools view  --threads $n_cores ~{prefix}.unaccounted_reads.bam > unaccounted_reads.bam
        if ~{default='true' reject_unaccounted} && [ -s unaccounted_reads.bam ]
            then exit 1; 
            else rm ~{prefix}.unaccounted_reads.bam
        fi 
        rm unaccounted_reads.bam
    >>>
 
    runtime {
        cpu: ncpu
        disk: disk_size + " GB"
        docker: 'stjudecloud/samtools:1.0.0'
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
        String outfilename = basename(bam, ".bam")+".flagstat.txt"
        Int max_retries = 1
        Int memory_gb = 5
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        samtools flagstat ~{bam} > ~{outfilename}
    >>>

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/samtools:1.0.0'
        memory: memory_gb + " GB"
        maxRetries: max_retries
    }

    output { 
       File outfile = outfilename
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool generates a FastQC quality control metrics report for the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate flagstat for"
    }
}

task index {
    input {
        File bam
        Int ncpu = 1
        String outfile = basename(bam)+".bai"
        Int max_retries = 1
        Int memory_gb = 15
        Boolean detect_nproc = false
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        n_cores=~{ncpu}
        if [ ${true='true' false='' detect_nproc} ]
        then
            n_cores=`nproc`
        fi

        samtools index -@ $n_cores ~{bam} ~{outfile}
    >>>

    runtime {
        cpu: ncpu
        disk: disk_size + " GB"
        docker: 'stjudecloud/samtools:1.0.0'
        memory: memory_gb + " GB"
        dx_instance_type: "azure:mem2_ssd1_x4"
        maxRetries: max_retries
    }

    output {
        File bai = outfile
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
        String outname = basename(bam, ".bam") + ".subsampled.bam"
        Int max_retries = 1
        Int desired_reads = 500000
        Int ncpu = 1
        Boolean detect_nproc = false
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ ${true='true' false='' detect_nproc} ]
        then
            n_cores=`nproc`
        fi

        if [[ "$(samtools view --threads $n_cores ~{bam} | head -n ~{desired_reads} | wc -l)" -ge "~{desired_reads}" ]]; then
            # the BAM has at least ~{desired_reads} reads, meaning we should
            # subsample it.
            initial_frac=0.00001
            initial_reads=$(samtools view --threads $n_cores -s $initial_frac ~{bam} | wc -l)
            frac=$( \
                awk -v desired_reads=~{desired_reads} \
                    -v initial_reads="$initial_reads" \
                    -v initial_frac=$initial_frac \
                        'BEGIN{printf "%1.8f", ( desired_reads / initial_reads * initial_frac )}' \
                )
            samtools view --threads $n_cores -h -b -s "$frac" ~{bam} > ~{outname}
        else
            # the BAM has less than ~{desired_reads} reads, meaning we should
            # just use it directly without subsampling.
            mv ~{bam} ~{outname}
        fi
    >>>

    output {
        File sampled_bam = outname
    }

    runtime {
        cpu: ncpu
        disk: disk_size + " GB"
        docker: 'stjudecloud/samtools:1.0.0'
        maxRetries: max_retries
    }
}
