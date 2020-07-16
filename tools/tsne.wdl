## # t-SNE
##
## This WDL tool implements a t-SNE (t-Distribute Stochastic Neighbor Embedding) plot to 
## visualize the relationship between input samples and a reference data set.
## (https://lvdmaaten.github.io/tsne/)

version 1.0

task plot {
    input {
        File counts
        Array[String]? inputs
        Array[File]? input_counts
        File blacklist
        File covariates
        File gencode_gtf
        String outfile = "tsne.html"
        String? tissue_type
        Int max_retries = 1
        Int memory_gb = 15
        Int disk_size = 15
    }

    command <<<
        mkdir counts
        tar --no-same-owner -C counts -zxf ${counts}
        count_arg=$(ls counts/*)

        itsne-main --debug-rscript --save-data -b ${blacklist} -g ${gencode_gtf} -c ${covariates} -o ${outfile} \
            ${true='--input-sample ' false='' defined(inputs)}${sep=' --input-sample ' inputs} \
            ${'--tissue-type ' + tissue_type} $count_arg ${sep=' ' input_counts}
    >>>
  

    runtime {
        disk: disk_size + " GB"
        memory: memory_gb + " GB"
        docker: 'stjudecloud/interactive-tsne:0.0.5'
        maxRetries: max_retries
    }


    output {
        File html = outfile
        File matrix = "tsne.txt"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
    }

    parameter_meta {
        covariates: "A tab delimited file with Sample, Diagnosis, Strandedness, LibraryType, ReadLength. Valid values for strandedness: [Stranded, Unstranded]. Valid values for LibraryType: [mRNA, totalRNA]. ReadLength should be specified similar to the form PE100bp."
    }
}

task append_input {
    input {
        Array[String] inputs
        File covariates_infile
        Int memory_gb = 5
    }

    command <<<
        set -euo pipefail
        
        cat ${covariates_infile} > "covariates.combined.tsv"
        for sample in ${sep=" " inputs}
        do
            echo -e "$sample\tinput\t$sample"
        done >> "covariates.combined.tsv"
    >>>

    runtime{
        memory: memory_gb + " GB"
    }

    output {
        File covariates_outfile = "covariates.combined.tsv"
    }

    meta {
        author: "Andrew Thrasher"
        email: "andrew.thrasher@stjude.org"
    }

    parameter_meta {
        inputs: "An array of input samples to append to the covariates file"
        covariates_infile: "A file of covariate data with the format 'sample\tdiagnosis\tstrandedness\tlibrary_type\tread_length'"
    }
}
