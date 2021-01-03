## # HTSeq
##
## This WDL tool wraps the [htseq](https://github.com/simon-anders/htseq) tool.
## HTSeq is a Python library for analyzing sequencing data.

version 1.0

task count {
    input {
        File bam
        File gtf
        String provided_strandedness
        String inferred_strandedness
        String outfile = basename(bam, ".bam") + ".feature-counts.txt"
        Int added_memory_gb = 20
        Int max_retries = 1
    }

    String stranded = if (provided_strandedness != "") then 
                        if (provided_strandedness == "Stranded-Reverse") then "reverse" else
                        if (provided_strandedness == "Stranded-Forward") then "yes" else
                        if (provided_strandedness == "Unstranded") then "no"
                        else "unknown-strand"
                      else 
                        if (inferred_strandedness == "Stranded-Reverse") then "reverse" else
                        if (inferred_strandedness == "Stranded-Forward") then "yes" else 
                        if (inferred_strandedness == "Unstranded") then "no" else
                        "unknown-strand" # this will intentionally cause htseq to error. You will need to manually specify
                                         # in this case

    Float bam_size = size(bam, "GiB")
    Float mem_size = bam_size + added_memory_gb
    Float gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(((bam_size + gtf_size) * 4) + 10)
 
    command {
        htseq-count -f bam \
            -r pos \
            -s ~{stranded} \
            -m union \
            -i gene_name \
            --secondary-alignments ignore \
            --supplementary-alignments ignore \
            ~{bam} \
            ~{gtf} \
            > ~{outfile}
    }

    runtime {
        memory: mem_size + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/htseq:1.0.0'
        maxRetries: max_retries
    }
   
    output {
        File out = "~{outfile}"
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool performs read counting for a set of features in the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
        gtf: "Input genomic features in GTF format to count reads for"
        added_memory_gb: "Amount of additional memory to add to the bam size"
    }
}

task calc_gene_lengths {
    input {
        File gtf
        String outfile = basename(gtf, ".gtf.gz") + ".genelengths.txt"
        Int max_retries = 1
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size = ceil(gtf_size * 2 + 10)

    command <<<
        GTF="~{gtf}" OUTFILE="~{outfile}" python - <<END
            import os
            import gtfparse

            gtf_name = os.environ['GTF']
            outfile = open(os.environ['OUTFILE'], 'w')

            gtf = gtfparse.parse_gtf(gtf_name)

            only_genes = gtf[gtf['feature'] == 'gene']
            only_exons = gtf[gtf['feature'] == 'exon']
            gene_start_offset = {}
            gene_end_offset = {}
            gene_exon_intersection = {}
            gene_total_exon_size = {}
            gene_length = {}

            for (index, value) in only_genes.iterrows():
                gene_name = value['gene_name']
                start = value['start']
                end = value['end']
                size = end - start
                
                if size <= 0:
                    raise RuntimeError("Size of gene is negative!")
                    
                gene_start_offset[gene_name] = start
                gene_end_offset[gene_name] = end
                gene_exon_intersection[gene_name] = np.zeros(size)
                gene_total_exon_size[gene_name] = 0
                gene_length[gene_name] = end - start
                
            for (index, value) in only_exons.iterrows():
                gene_name = value['gene_name']
                offset = gene_start_offset[gene_name]
                start = value['start'] - offset
                end = value['end'] - offset
                exon_length = end - start
                gene_exon_intersection[gene_name][start:end] = 1
                gene_total_exon_size[gene_name] += exon_length

            results = []
            print("Gene name\tlength", file=outfile)
            for (gene, exonic_intersection) in gene_exon_intersection.items():
                length = np.sum(exonic_intersection).astype(int)
                print(f"{gene}\t{length}")
        END
    >>>

    runtime {
        memory: "8 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/gtfparse:branch-ESTIMATE-1.0.0'
        maxRetries: max_retries
    }

    output {
        File out = "~{outfile}"
    }
}