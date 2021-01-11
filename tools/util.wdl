## # Utilities
##
## This WDL tool includes custom scripts to parse and reformat 
## task output as part of a workflow. 

version 1.0

task get_read_groups {
    input {
        File bam
        Int max_retries = 1
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size = ceil((bam_size * 2) + 10)

    command <<<
        set -euo pipefail
        
        samtools view -H ~{bam} | grep "@RG" \
            | cut -f 2- \
            | sed -e 's/\t/ /g' \
            | awk '{print}' ORS=' , ' \
            | sed 's/ , $//' >> stdout.txt
    >>>

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/samtools:1.0.0'
        maxRetries: max_retries
    }

    output { 
        String out = read_string("stdout.txt")
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool is a utility to get read group information from a BAM file and write it out to as a string" 
    }

    parameter_meta {
        bam: "Input BAM format file to get read groups from"
    }
}

task file_prefix {
    input {
        File in_file
        Int max_retries = 1
    }

    Float infile_size = size(in_file, "GiB")
    Int disk_size = ceil(infile_size * 2)

    command <<<
        set -euo pipefail
        
        basename ~{in_file} | awk -F '.' '{print $1}' > stdout.txt
    >>>

    runtime {
        disk: disk_size + " GB"
        docker: 'stjudecloud/util:1.0.0'
        maxRetries: max_retries
    }

    output { 
        String out = read_string("stdout.txt")
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
import os  # lint-check: ignore
import gtfparse  # lint-check: ignore
import numpy as np  # lint-check: ignore

gtf_name = os.environ['GTF']
outfile = open(os.environ['OUTFILE'], 'w')

gtf = gtfparse.read_gtf(gtf_name)

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
for (gene, exonic_intersection) in sorted(gene_exon_intersection.items()):
    length = np.sum(exonic_intersection).astype(int)
    print(f"{gene}\t{length}", file=outfile)
END
    >>>

    runtime {
        memory: "8 GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/gtfparse:1.0.0'
        maxRetries: max_retries
    }

    output {
        File out = "~{outfile}"
    }
}
