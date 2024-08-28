version 1.1

task pre {
    meta {
        description: "Run Juicer pre to generate the .hic file"
        outputs: {
            hic: "Reads in .hic format"
        }
    }

    parameter_meta {
        pairs: "Read pairs file"
        genome_id: {
            description: "Genome ID",
            choices: [
                "hg18",
                "hg19",
                "hg38",
                "dMel",
                "mm9",
                "mm10",
                "anasPlat1",
                "bTaurus3",
                "canFam3",
                "equCab2",
                "galGal4",
                "Pf3D7",
                "sacCer3",
                "sCerS288c",
                "susScr3",
                "TAIR10",
            ],
        }
        prefix: "Prefix for output files. The extension `.hic` will be added."
        diagonal: "Only calculate intra chromosome (diagonal)"
        min_count: "Only write cells with count above threshold"
        restriction_sites: {
            description: "Calculate fragment map. Requires restriction site file; each line should start with the chromosome name followed by the position of each restriction site on that chromosome, in numeric order, and ending with the size of the chromosome.",
            external_help: "https://github.com/aidenlab/juicer/wiki/Pre#restriction-site-file-format",
        }
        stats_file: {
            description: "Add the text statistics file to the Hi-C file header",
            external_help: "https://github.com/theaidenlab/juicer/blob/master/UGER/scripts/statistics.pl",
        }
        graphs_file: {
            description: "Add the text graphs file to the Hi-C file header",
            external_help: "https://github.com/theaidenlab/juicer/blob/master/UGER/scripts/statistics.pl",
        }
        resolutions: "Only calculate specific resolutions"
        maq_filter: "Filter by MAPQ score greater than or equal to this value"
        chromosome_filter: "Only calculate map on specific chromosome"
        disable_normalize: "Don't normalize the matrices"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File pairs
        File? restriction_sites
        File? stats_file
        File? graphs_file
        Array[Int]? resolutions
        String? chromosome_filter
        Int? maq_filter
        String genome_id = "hg38"
        String prefix = basename(pairs, ".bsorted.pairs.gz")
        Boolean diagonal = false
        Boolean disable_normalize = false
        Int min_count = 0
        Int modify_disk_size_gb = 0
    }

    Int disk_size_gb = ceil(size(pairs, "GiB") * 2 ) + modify_disk_size_gb

    command <<<
        juicer_tools \
            pre \
            ~{if diagonal then "-d" else ""} \
            ~{if defined(restriction_sites) then "-f " + restriction_sites else ""} \
            -m ~{min_count} \
            ~{if defined(maq_filter) then "-q " + maq_filter else ""} \
            ~{if defined(chromosome_filter) then "-c " + chromosome_filter else ""} \
            ~{(
                if defined(resolutions)
                then "-r " + sep(",", select_first([resolutions, []]))
                else ""
            )} \
            ~{if defined(stats_file) then "-s " + stats_file else ""} \
            ~{if defined(graphs_file) then "-g " + graphs_file else ""} \
            ~{if disable_normalize then "-n" else ""} \
            ~{pairs} \
            ~{prefix}.hic \
            ~{genome_id}
    >>>

    output {
        File hic = "~{prefix}.hic"
    }

    runtime {
        cpu: 1
        memory: "14 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/juicer:branch-hic_workflow-1.0.13-0"
        maxRetries: 1
    }
}
