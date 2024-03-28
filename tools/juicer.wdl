version 1.1

task pre {
    meta {
        description: "Run Juicer pre to generate the .hic file"
        output: {
            hic: "Reads in .hic format"
        }
    }

    parameter_meta {
        pairs: "Read pairs file"
        genomeID: {
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
                "TAIR10"
            ],
        }
        prefix: "Prefix for output files. The extension `.hic` will be added."
        diagonal: "Only calculate intra chromosome (diagonal)"
        min_count: "Only write cells with count above threshold"
        restriction_sites: {
            description: "Calculate fragment map. Requires restriction site file; each line should start with the chromosome name followed by the position of each restriction site on that chromosome, in numeric order, and ending with the size of the chromosome."
            external_help: "https://github.com/aidenlab/juicer/wiki/Pre#restriction-site-file-format"
        }
        stats_file: {
            description: "Add the text statistics file to the Hi-C file header"
            external_help: "https://github.com/theaidenlab/juicer/blob/master/UGER/scripts/statistics.pl"
        }
        graphs_file: {
            description: "Add the text graphs file to the Hi-C file header"
            external_help: "https://github.com/theaidenlab/juicer/blob/master/UGER/scripts/statistics.pl"
        }
        resolutions: "Only calculate specific resolutions"
        maq_filter: "Filter by MAPQ score greater than or equal to this value"
        chromosome_filter: "Only calculate map on specific chromosome"
        disable_normalize: "Don't normalize the matrices"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File pairs
        String genomeID = "hg38"
        String prefix = basename(pairs, ".pairs.txt")
        Boolean diagonal = false
        Int min_count = 0
        Int modify_disk_size_gb = 0
        File? restriction_sites
        File? stats_file
        File? graphs_file
        Array[Int]? resolutions
        Int? maq_filter
        String? chromosome_filter
        Boolean disable_normalize = false
    }

    Int disk_size_gb = ceil(size(pairs, "GiB") * 2 ) + modify_disk_size_gb 

    command <<<
        juicer_tools \
            pre \
            -d ~{diagonal} \
            ~{if defined(restriction_sites) then "-f " + restriction_sites else ""} \
            -m ~{min_count} \
            ~{if defined(maq_filter) then "-q " + maq_filter else ""} \
            ~{if defined(chromosome_filter) then "-c " + chromosome_filter else ""} \
            ~{if defined(resolutions) then "-r " + sep(',', select_first([resolutions, []])) else ""} \
            ~{if defined(stats_file) then "-s " + stats_file else ""} \
            ~{if defined(graphs_file) then "-g " + graphs_file else ""} \
            ~{if disable_normalize then "-n" else ""} \
            ~{pairs} \
            ~{prefix}.hic \
            ~{genomeID}
    >>>

    output {
        File hic = "~{prefix}.hic"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: "aidenlab/juicer:1.0.13"
        maxRetries: 1
    }
}