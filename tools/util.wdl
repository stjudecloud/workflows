## # Utilities

version 1.1

task download {
    meta {
        description: "Uses wget to download a file from a remote URL to the local filesystem"
        outputs: {
            downloaded_file: "File downloaded from provided URL"
        }
    }

    parameter_meta {
        url: "URL of the file to download"
        outfile_name: "Name of the output file"
        disk_size_gb: "Disk space to allocate for task, specified in GB"
        md5sum: "Optional md5sum to check against downloaded file. Recommended to use in order to catch corruption or an unintentional file swap."
    }

    input {
        String url
        String outfile_name
        Int disk_size_gb
        String? md5sum
    }

    command <<<
        set -euo pipefail

        wget ~{url} -O ~{outfile_name}

        if [ -n "~{md5sum}" ]; then
            echo "~{md5sum}  ~{outfile_name}" > ~{outfile_name}.md5
            md5sum -c ~{outfile_name}.md5
        fi
    >>>

    output {
        File downloaded_file = outfile_name
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:branch-scripts-2.1.0"
        maxRetries: 1
    }
}

task get_read_groups {
    meta {
        description: "Gets read group information from a BAM file and writes it out to as a string"
        outputs: {
            read_groups: "An array of strings containing read group information. If `clean = true`, the `@RG\t` prefix is stripped and tabs are replaced with spaces. If `clean = false`, each unmodified @RG line will be its own entry in output array `read_groups`."
        }
    }

    parameter_meta {
        bam: {
            description: "Input BAM format file to get read groups from",
            stream: true,
        }
        clean: {
            description: "Clean @RG lines to remove the `@RG\t` prefix and use spaces instead of tabs (true) or output @RG lines of the header without further processing (false)?",
            help: "`clean = true` output matches the formatting of the `read_group_to_string` task in `../data_structures/read_group.wdl`",
            group: "common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        Boolean clean = true
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        if ~{clean}; then
            samtools view -H ~{bam} \
                | grep "^@RG" \
                | cut -f 2- \
                | sed -e 's/\t/ /g' \
                > read_groups.txt
        else
            samtools view -H ~{bam} | grep "^@RG" > read_groups.txt
        fi
    >>>

    output {
        Array[String] read_groups = read_lines("read_groups.txt")
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task split_string {
    # Currently (v1.1) no way to do this using the WDL standard library.
    # Revisit task in future version updates, can hopefully be replaced.
    meta {
        description: "Split a string into an array of strings based on a delimiter"
        outputs: {
            split_strings: "Split string as an array"
        }
    }

    parameter_meta {
        string: "String to split on occurences of `delimiter`"
        delimiter: {
            description: "Delimiter on which to split `input_string`",
            group: "common",
        }
    }

    input {
        String string
        String delimiter = " , "
    }

    command <<<
        set -euo pipefail

        echo ~{string} | sed 's/~{delimiter}/\n/g' > split_strings.txt
    >>>

    output {
        Array[String] split_strings = read_lines("split_strings.txt")
    }

    runtime {
        memory: "4 GB"
        disks: "10 GB"
        container: "ghcr.io/stjudecloud/util:branch-scripts-2.1.0"
        maxRetries: 1
    }
}

task calc_gene_lengths {
    meta {
        description: "Calculate gene lengths from a GTF feature file using the non-overlapping exonic length algorithm"
        help: "The non-overlapping exonic length algorithm can be implemented as the sum of each base covered by at least one exon; where each base is given a value of 1 regardless of how many exons overlap it."
        outputs: {
            gene_lengths: "A two column headered TSV file with gene names in the first column and feature lengths (as integers) in the second column"
        }
    }

    parameter_meta {
        gtf: "GTF feature file"
        outfile_name: "Name of the gene lengths file"
        idattr: {
            description: "GTF attribute to be used as feature ID. The value of this attribute will be used as the first column in the output file.",
            group: "common",
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File gtf
        String outfile_name = basename(gtf, ".gtf.gz") + ".genelengths.txt"
        String idattr = "gene_name"
        Int modify_disk_size_gb = 0
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size_gb = ceil(gtf_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        GTF="~{gtf}" OUTFILE="~{outfile_name}" IDATTR="~{idattr}" python - <<END
        import os  # lint-check: ignore
        import gtfparse  # lint-check: ignore
        import numpy as np  # lint-check: ignore
        from collections import defaultdict  # lint-check: ignore

        gtf_name = os.environ["GTF"]
        outfile = open(os.environ["OUTFILE"], "w")
        id_attr = os.environ["IDATTR"]

        gtf = gtfparse.read_gtf(gtf_name)

        only_exons = gtf[gtf["feature"] == "exon"]
        exon_starts = defaultdict(lambda: [])
        exon_ends = defaultdict(lambda: [])
        gene_start_offset = {}
        gene_end_offset = {}
        gene_exon_intersection = {}

        for _index, value in only_exons.iterrows():
            feature_id = value[id_attr]
            start = value["start"]
            end = value["end"] + 1  # end is inclusive in GTF
            exon_starts[feature_id].append(start)
            exon_ends[feature_id].append(end)
            if feature_id not in gene_start_offset:
                gene_start_offset[feature_id] = start
                gene_end_offset[feature_id] = end
            else:
                gene_start_offset[feature_id] = min(gene_start_offset[feature_id], start)
                gene_end_offset[feature_id] = max(gene_end_offset[feature_id], end)

        for feature_id in exon_starts:
            gene_exon_intersection[feature_id] = np.full(
                gene_end_offset[feature_id] - gene_start_offset[feature_id], False
            )

            for start, end in zip(exon_starts[feature_id], exon_ends[feature_id]):
                gene_exon_intersection[feature_id][
                    start - gene_start_offset[feature_id]
                    : end - gene_start_offset[feature_id]
                ] = True

        print("feature\tlength", file=outfile)
        for gene, exonic_intersection in sorted(gene_exon_intersection.items()):
            # np.count_nonzero() is faster than sum
            # np.count_nonzero() evaluates the "truthfulness" of
            # of all elements (by calling their '.__bool__()' method)
            length = np.count_nonzero(exonic_intersection)
            print(f"{gene}\t{length}", file=outfile)

        outfile.close()

        END
    >>>

    output {
        File gene_lengths = "~{outfile_name}"
    }

    runtime {
        memory: "16 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0"
        maxRetries: 1
    }
}

task compression_integrity {
    meta {
        description: "Checks the compression integrity of a bgzipped file"
        outputs: {
            check: "Dummy output to indicate success and to enable call-caching"
        }
    }

    parameter_meta {
        bgzipped_file: "Input bgzipped file to check integrity of"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bgzipped_file
        Int modify_disk_size_gb = 0
    }

    Float file_size = size(bgzipped_file, "GiB")
    Int disk_size_gb = ceil(file_size) + 10 + modify_disk_size_gb

    command <<<
        bgzip -t ~{bgzipped_file}
    >>>

    output {
        String check = "passed"
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task add_to_bam_header {
    meta {
        description: "Adds another line of text to the bottom of a BAM header"
        outputs: {
            reheadered_bam: "The BAM after its header has been modified"
        }
    }

    parameter_meta {
        bam: "Input BAM format file which will have its header added to"
        additional_header: "A string to add as a new line in the BAM header. No format checking is done, so please ensure you do not invalidate your BAM with this task. Add only spec compliant entries to the header."
        prefix: "Prefix for the reheadered BAM. The extension `.bam` will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String additional_header
        String prefix = basename(bam, ".bam") + ".reheader"
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".bam"

    command <<<
        samtools view -H ~{bam} > header.sam
        echo "~{additional_header}" >> header.sam
        samtools reheader -P header.sam ~{bam} > ~{outfile_name}
    >>>

    output {
        File reheadered_bam = outfile_name
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
        maxRetries: 1
    }
}

task unpack_tarball {
    meta {
        description: "Accepts a `.tar.gz` archive and converts it into a flat array of files. Any directory structure of the archive is ignored."
        outputs: {
            tarball_contents: "An array of files found in the input tarball"
        }
    }

    parameter_meta {
        tarball: "A `.tar.gz` archive to unpack into individual files"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File tarball
        Int modify_disk_size_gb = 0
    }

    Float tarball_size = size(tarball, "GiB")
    Int disk_size_gb = ceil(tarball_size * 8) + modify_disk_size_gb

    command <<<
        set -euo pipefail

        mkdir unpacked_tarball
        tar -C unpacked_tarball -xzf ~{tarball} --no-same-owner
        # pipe through sort because otherwise order is random (dependent on filesystem)
        find unpacked_tarball/ -type f | LC_ALL=C sort > file_list.txt
    >>>

    output {
        Array[File] tarball_contents = read_lines("file_list.txt")
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:branch-scripts-2.1.0"
        maxRetries: 1
    }
}

task make_coverage_regions_bed {
    meta {
        description: "Takes in a GTF file, converts it to BED, then filters it down to a 3 column BED file from all lines which match a given feature type"
        outputs: {
            bed: "3 column BED file corresponding to all records in the input GTF with a feature type matching `feature_type`",
        }
    }

    parameter_meta {
        gtf: "Gzipped GTF feature file from which to derive a coverage regions BED file"
        feature_type: {
            description: "Feature type to filter on. Only lines with this feature type will be included in the output BED file.",
            help: "`choices` below are the possible values from a GENCODE GTF file. If you are using a different GTF source, you may need to adjust this parameter.",
            choices: [
                "gene",
                "transcript",
                "exon",
                "CDS",
                "UTR",
                "start_codon",
                "stop_codon",
                "Selenocysteine"
            ],
        }
        outfile_name: "Name of the output BED file"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File gtf
        String feature_type
        String outfile_name = basename(gtf, "gtf.gz") + feature_type + ".bed"
        Int modify_disk_size_gb = 0
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size_gb = ceil(gtf_size * 1.2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        gunzip -c ~{gtf} \
            | gtf2bed \
            | awk '$8 == "~{feature_type}" {print $1 "\t" $2 "\t" $3}' \
            > ~{outfile_name}
    >>>

    output {
        File bed = outfile_name
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0"
        maxRetries: 1
    }
}

task global_phred_scores {
    meta {
        description: "Calculates statistics about PHRED scores of the input BAM"
        outputs: {
            phred_scores: "Headered TSV file containing PHRED score statistics"
        }
    }

    parameter_meta {
        bam: "Input BAM format file to calculate PHRED score statistics for"
        prefix: "Prefix for the output TSV file. The extension `.global_PHRED_scores.tsv` will be added."
        fast_mode: "Enable fast mode (true) or calculate statistics for *_every_* base in the BAM (false)?"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        String prefix = basename(bam, ".bam")
        Boolean fast_mode = true
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    String outfile_name = prefix + ".global_PHRED_scores.tsv"

    #@ except: LineWidth
    command <<<
        set -euo pipefail

        python3 /scripts/util/calc_global_phred_scores.py \
            ~{if fast_mode then "--fast_mode" else ""} \
            ~{bam} \
            ~{prefix}
    >>>

    output {
        File phred_scores = "~{outfile_name}"
    }

    runtime {
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:branch-scripts-2.1.0"
        maxRetries: 1
    }
}

task qc_summary {
    meta {
        description: "**[OUT OF DATE]** This WDL task pulls out keys metrics that can provide a high level overview of the sample, without needing to examine the entire MultiQC report. Currently, these key metrics come from Qualimap and ngsderive."
        outputs: {
            summary: "QC summary file in JSON format"
        }
    }

    parameter_meta {
        multiqc_tar_gz: "MultiQC report tarball from which to extract key metrics"
        outfile_name: "Name for the JSON file"
    }

    input {
        File multiqc_tar_gz
        String outfile_name
            = basename(multiqc_tar_gz, ".multiqc.tar.gz") + ".qc_summary.json"
    }

    String sample_name = basename(multiqc_tar_gz, ".multiqc.tar.gz")

    #@ except: LineWidth
    command <<<
        set -euo pipefail

        tar -xzf "~{multiqc_tar_gz}" --no-same-owner
        gen_stats_file=~{sample_name}.multiqc/multiqc_data/multiqc_general_stats.txt

        TOTAL_READS=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-total_reads $gen_stats_file | tail -n 1 | awk '{ printf("%.0f", $1) }')
        PERCENT_ALIGNED=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-percentage_aligned $gen_stats_file | tail -n 1 | awk '{ printf("%.3f", $1) }')
        MEAN_COVERAGE=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-mean_coverage $gen_stats_file | tail -n 1 | awk '{ printf("%.3f", $1) }')
        INSERT_SIZE=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-median_insert_size $gen_stats_file | tail -n 1)
        MEAN_GC_CONTENT=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-avg_gc $gen_stats_file | tail -n 1 | awk '{ printf("%.3f", $1) }')
        READ_LENGTH=$(csvcut -t -c ngsderive_mqc-generalstats-ngsderive-consensusreadlength $gen_stats_file | tail -n 1)
        PLATFORM=$(csvcut -t -c ngsderive_mqc-generalstats-ngsderive-instrument $gen_stats_file | tail -n 1)

        THIRTYX_PERCENT=$(csvcut -t -c QualiMap_mqc-generalstats-qualimap-30_x_pc $gen_stats_file | tail -n 1)
        DUP_PERCENT=$(csvcut -t -c FastQC_mqc-generalstats-fastqc-percent_duplicates $gen_stats_file | tail -n 1)

        STRANDEDNESS=$({ csvcut -t -c ngsderive_mqc-generalstats-ngsderive-predicted $gen_stats_file || echo "Not Applicable" ; } | tail -n 1)

        jq -n \
            --arg TOTAL_READS "$TOTAL_READS" \
            --arg PERCENT_ALIGNED "$PERCENT_ALIGNED" \
            --arg MEAN_COVERAGE "$MEAN_COVERAGE" \
            --arg INSERT_SIZE "$INSERT_SIZE" \
            --arg MEAN_GC_CONTENT "$MEAN_GC_CONTENT" \
            --arg READ_LENGTH "$READ_LENGTH" \
            --arg PLATFORM "$PLATFORM" \
            --arg THIRTYX_PERCENT "$THIRTYX_PERCENT" \
            --arg DUP_PERCENT "$DUP_PERCENT" \
            --arg STRANDEDNESS "$STRANDEDNESS" \
            '{
                total_reads: ($TOTAL_READS | tonumber),
                percent_aligned: ($PERCENT_ALIGNED | tonumber),
                mean_coverage: ($MEAN_COVERAGE | tonumber),
                median_insert_size: ($INSERT_SIZE | tonumber),
                mean_percent_GC: ($MEAN_GC_CONTENT | tonumber),
                inferred_read_length: ($READ_LENGTH | tonumber),
                inferred_sequencing_platform: $PLATFORM,
                percent_thirtyX_coverage: ($THIRTYX_PERCENT | tonumber),
                percent_duplicate: ($DUP_PERCENT | tonumber),
                inferred_strandedness: $STRANDEDNESS
            }' > ~{outfile_name}
    >>>

    output {
        File summary = "~{outfile_name}"
    }

    runtime {
        memory: "4 GB"
        disks: "10 GB"
        container: "ghcr.io/stjudecloud/util:branch-scripts-2.1.0"
        maxRetries: 1
    }
}

task split_fastq {
    meta {
        description: "Splits a FASTQ into multiple files based on the number of reads per file"
        outputs: {
            fastqs: "Array of FASTQ files, each containing a subset of the input FASTQ"
        }
    }

    parameter_meta {
        fastq: {
            description: "Gzipped FASTQ file to split",
            stream: true,
        }
        reads_per_file: "Number of reads to include in each output FASTQ file"
        prefix: "Prefix for the FASTQ files. The extension `.fastq.gz` (preceded by a split index) will be added."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
        ncpu: "Number of cores to allocate for task"
    }

    input {
        File fastq
        String prefix = sub(
            basename(fastq),
            "(fastq|fq)\\.gz$",
            ""
        )
        Int reads_per_file = 10000000
        Int modify_disk_size_gb = 0
        Int ncpu = 2
    }

    Float fastq_size = size(fastq, "GiB")
    Int disk_size_gb = ceil(fastq_size * 5) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        let "lines = ~{reads_per_file} * 4"
        zcat ~{fastq} | split -l $lines -d -a 6 - ~{prefix}

        for file in "~{prefix}"*; do
            mv "$file" "${file}.fastq"
            echo "gzip ${file}.fastq" >> cmds
        done

        parallel --jobs ~{ncpu} < cmds
    >>>

    output {
        Array[File] fastqs = glob("~{prefix}*")
    }

    runtime {
        cpu: ncpu
        memory: "4 GB"
        disks: "~{disk_size_gb} GB"
        container: "ghcr.io/stjudecloud/util:branch-scripts-2.1.0"
        maxRetries: 1
    }
}
