version 1.1

task trim_single_end {
    meta {
        description: "Performs single end trimming of reads using the Trimmomatic tool."
        outputs: {
            trimmed_reads: "The trimmed single-end reads in FASTQ format."
        }
    }

    parameter_meta {
        input_file: "Input FASTQ format file to run Trimmomatic on."
        output_name: "Output FASTQ format file of trimmed reads."
        illumina_clip: "Cut adapter and other illumina-specific sequences from the read."
        leading: "Cut bases off the start of a read, if below a threshold quality"
        threads: "Number of threads to use for trimming."
        trailing: "Cut bases off the end of a read, if below a threshold quality."
        window_size: "Specifies the number of bases to average across for the sliding window."
        window_quality: "Specifies the average quality required in the sliding window."
        minlen: "Drop the read if it is below a specified length."
        crop: "Cut the read to a specified length."
        headcrop: "Cut the specified number of bases from the start of the read"
        phred64: "Convert quality scores to Phred-64"
    }

    input {
        File input_file
        String outfile_name = (
            basename(input_file, ".fastq.gz") + ".trimmed.fastq.gz"
        )
         String prefix = sub(
            basename(input_file),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        String illumina_clip
        Int? crop
        Int? headcrop
        Int? leading
        Int minlen
        Int threads = 8
        Int? trailing
        Int window_size
        Int window_quality
        Boolean phred64
    }

    command <<<
        ls -l /usr/local/share/trimmomatic-0.36-5/adapters/
        trimmomatic SE -threads ~{threads} \
        -trimlog trimlog.txt ~{input_file} output.fq.gz \
        ~{if (phred64) then "-phred64" else "-phred33"} \
        ILLUMINACLIP:~{illumina_clip} \
        ~{if defined(leading) then "LEADING:~{leading}" else ""} \
        ~{if defined(trailing) then "TRAILING:~{trailing}" else ""} \
        SLIDINGWINDOW:~{window_size}:~{window_quality} \
        MINLEN:~{minlen} \
        ~{if defined(crop) then "CROP:~{crop}" else ""} \
        ~{if defined(headcrop) then "HEADCROP:~{headcrop}" else ""} \
        ; mv output.fq.gz ~{outfile_name}
    >>>

    output {
        File trimmed_reads = "~{outfile_name}"
    }

    runtime {
        memory: "16 GB"
        disks: "100 GB"
        container: "quay.io/biocontainers/trimmomatic:0.36--5"
        cpu: 8
        maxRetries: 0
    }
}

task trim_paired_end {
    meta {
        description: "Performs paired end trimming of reads using the Trimmomatic tool."
        outputs: {
            trimmed_reads_forward: "The forward trimmed paired-end read in FASTQ format."
            trimmed_reads_reverse: "The reverse trimmed paired-end read in FASTQ format."
        }
    }

    parameter_meta {
        input_file_forward: "Forward input FASTQ format file to run Trimmomatic on."
        input_file_reverse: "Reverse input FASTQ format file to run Trimmomatic on."
        output_name: "Output FASTQ format file of trimmed reads."
        leading: "Cut bases off the start of a read, if below a threshold quality."
        trailing: "Cut bases off the end of a read, if below a threshold quality."
        window_size: "Specifies the number of bases to average across for the sliding window."
        window_quality: "Specifies the average quality required in the sliding window."
        minlen: "Drop the read if it is below a specified length."
    }

    input {
        File input_file_forward
        File input_file_reverse
        String outfile_name_forward = (
            basename(input_file_forward, ".fastq.gz") + ".forward.trimmed.fastq.gz"
        )
        String outfile_name_reverse = (
            basename(input_file_reverse, ".fastq.gz") + ".reverse.trimmed.fastq.gz"
        )
         String prefix = sub(
            basename(input_file_forward),
            "([_\\.][rR][12])?(\\.subsampled)?\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        String illumina_clip
        Int? crop
        Int? headcrop
        Int? leading
        Int minlen
        Int threads = 8
        Int? trailing
        Int window_size
        Int window_quality
        Boolean phred64
    }

    command <<<
        ls -l /usr/local/share/trimmomatic-0.36-5/adapters/
        trimmomatic PE -phred33 -threads ~{threads} \
        -trimlog trimlog.txt ~{input_file_forward} ~{input_file_reverse} \
        output_forward.fq.gz output_reverse.fq.gz \
        ~{if (phred64) then "-phred64" else "-phred33"} \
        ILLUMINACLIP:~{illumina_clip} \
        ~{if defined(leading) then "LEADING:~{leading}" else ""} \
        ~{if defined(trailing) then "TRAILING:~{trailing}" else ""} \
        SLIDINGWINDOW:~{window_size}:~{window_quality} \
        MINLEN:~{minlen} \
        ~{if defined(crop) then "CROP:~{crop}" else ""} \
        ~{if defined(headcrop) then "HEADCROP:~{headcrop}" else ""} \
        ; mv output_forward.fq.gz ~{outfile_name_forward} \
        ; mv output_reverse.fq.gz ~{outfile_name_reverse}
    >>>

    output {
        File forward_trimmed_reads = "~{outfile_name_forward}"
        File reverse_trimmed_reads = "~{outfile_name_reverse}"
    }

    runtime {
        memory: "16 GB"
        disks: "100 GB"
        container: "quay.io/biocontainers/trimmomatic:0.36--5"
        cpu: 8
        maxRetries: 0
    }
}