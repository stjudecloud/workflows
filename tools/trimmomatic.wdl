version 1.1

task trim {
    meta {
        description: "Performs paired end trimming of reads using the Trimmomatic tool."
        outputs: {
            trimmed_read_one: "The forward trimmed paired-end read in FASTQ format."
            trimmed_read_two: "The reverse trimmed paired-end read in FASTQ format."
        }
    }

    parameter_meta {
        read_one: "Forward input FASTQ format file to run Trimmomatic on."
        read_two: "Reverse input FASTQ format file to run Trimmomatic on."
        output_name: "Output FASTQ format file of trimmed reads."
        leading: "Cut bases off the start of a read, if below a threshold quality."
        trailing: "Cut bases off the end of a read, if below a threshold quality."
        window_size: "Specifies the number of bases to average across for the sliding window."
        window_quality: "Specifies the average quality required in the sliding window."
        minlen: "Drop the read if it is below a specified length."
    }

    input {
        File read_one
        File? read_two
        String? illumina_clip
        String prefix = sub(
            basename(read_one),
            "([_\\.][rR][12])?[_\\.](.*)\\.(fastq|fq)(\\.gz)?$",
            ""
        )
        String outfile_name_one = "~{prefix}.R1.trimmed.fastq.gz"
        String outfile_name_two = "~{prefix}.R2.trimmed.fastq.gz"
        Int? crop
        Int? headcrop
        Int? leading
        Int minlen
        Int threads = 8
        Int? trailing
        Int window_quality
        Int window_size
        Boolean phred64 = false
    }

    command <<<
        ls -l /usr/local/share/trimmomatic-0.36-5/adapters/
        trimmomatic ~{if defined(read_two) then "PE" else "SE"} -threads ~{threads} \
        -trimlog trimlog.txt ~{read_one} ~{if defined(read_two) then "~{read_two}" else ""} \
        ~{outfile_name_one} ~{if defined(read_two) then "~{outfile_name_two}" else ""} \
        ~{if defined(phred64) && phred64 then "-phred64" else ""} \
        ~{if defined(illumina_clip) then "ILLUMINACLIP:~{illumina_clip}" else ""} \
        ~{if defined(leading) then "LEADING:~{leading}" else ""} \
        ~{if defined(trailing) then "TRAILING:~{trailing}" else ""} \
        SLIDINGWINDOW:~{window_size}:~{window_quality} \
        MINLEN:~{minlen} \
        ~{if defined(crop) then "CROP:~{crop}" else ""} \
        ~{if defined(headcrop) then "HEADCROP:~{headcrop}" else ""} 
    >>>

    output {
        File trimmed_read_one = "~{outfile_name_one}"
        File? trimmed_read_two = "~{outfile_name_two}"
    }

    runtime {
        memory: "16 GB"
        disks: "100 GB"
        container: "quay.io/biocontainers/trimmomatic:0.36--5"
        cpu: 8
        maxRetries: 0
    }
}