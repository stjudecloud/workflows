version 1.1

task trim_single_end {
    meta {
        description: "Performs both single end and paired end trimming of reads using the Trimmomatic tool."
        outputs: {
            trimmed_reads: "The trimmed reads in FASTQ format."
        }
    }

    parameter_meta {
        input_file: "Input FASTQ format file to run Trimmomatic on."
        output_name: "Output FASTQ format file of trimmed reads."
        # adapter_file: "Adapter file to use for trimming."
        leading: "Cut bases off the start of a read, if below a threshold quality."
        trailing: "Cut bases off the end of a read, if below a threshold quality."
        window_size: "Specifies the number of bases to average across for the sliding window."
        window_quality: "Specifies the average quality required in the sliding window."
        minlen: "Drop the read if it is below a specified length."
    }

    input {
        File inputfile
        String output_name
        # String adapter_file
        Int leading
        Int trailing
        Int window_size
        Int window_quality
        Int minlen
    }

    command <<<
        trimmomatic SE -phred33 -threads 8 \
        -trimlog trimlog.txt ~{inputfile} output.fq.gz \
        ILLUMINACLIP:/usr/local/share/trimmomatic-0.36-5/adapters/TruSeq3-SE.fa:2:30:10 \
        LEADING:~{leading} \
        TRAILING:~{trailing} \
        SLIDINGWINDOW:~{window_size}:~{window_quality} \
        MINLEN:~{minlen} 
        ls -l
    >>>

    output {
        File trimmed_reads = "output.fq.gz"
    }

    runtime {
        memory: "16 GB"
        disks: "100 GB"
        container: "quay.io/biocontainers/trimmomatic:0.36--5"
        cpu: 8
        maxRetries: 0
    }
}