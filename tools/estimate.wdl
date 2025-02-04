## [Homepage](https://bioinformatics.mdanderson.org/estimate/)

version 1.1

task run_estimate {
    meta {
        description: "**[DEPRECATED]** Given a gene expression file, run the ESTIMATE software package"
        outputs: {
            estimate_file: "The results file of the ESTIMATE software package"
        }
        deprecated: true
    }

    parameter_meta {
        gene_expression_file: "A 2 column headered TSV file with 'Gene name' in the first column and gene expression values (as floats) in the second column. Can be generated with the `calc_tpm` task."
        outfile_name: "Name of the ESTIMATE output file"
        memory_gb: "RAM to allocate for task, specified in GB"
        disk_size_gb: "Disk space to allocate for task, specified in GB"
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        File gene_expression_file
        String outfile_name = (
            basename(gene_expression_file, ".TPM.txt") + ".ESTIMATE.gct"
        )
        Int memory_gb = 4
        Int disk_size_gb = 10
        Int max_retries = 1
    }

    #@ except: LineWidth
    command <<<
        cp "~{gene_expression_file}" gene_expression.txt
        Rscript - <<END
        library("estimate")

        infile <- read.table(file = "gene_expression.txt", sep = '\t', header = TRUE)
        filtered <- infile[infile$"Gene.name" %in% common_genes[['GeneSymbol']], ]
        write.table(filtered, sep = "\t", file = "filtered.tsv", row.names = FALSE, quote = FALSE)
        outputGCT("filtered.tsv", "gene_expression.gct")
        estimateScore("gene_expression.gct", "common_estimate.gct", platform = "illumina")
        END
        mv common_estimate.gct "~{outfile_name}"
    >>>

    output {
        File estimate_file = "~{outfile_name}"
    }

    runtime {
        memory: "~{memory_gb} GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/estimate:2.0.0"
        maxRetries: max_retries
    }
}
