## # Utilities
#
# SPDX-License-Identifier: MIT
# Copyright St. Jude Children's Research Hospital
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
        disk: "~{disk_size_gb} GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}

task get_read_groups {
    meta {
        description: "Gets read group information from a BAM file and writes it out to as a string"
        outputs: {
            read_groups: "An array of strings containing read group information. If `format_for_star = true`, all found read groups are contained in one string (`read_groups[0]`). If `format_for_star = false`, each found @RG line will be its own entry in output array `read_groups`."
        }
    }

    parameter_meta {
        bam: "Input BAM format file to get read groups from"
        format_for_star: {
            description: "Format read group information for the STAR aligner (true) or output @RG lines of the header without further processing (false)? STAR formatted results will be an array of length 1, where all found read groups are contained in one string (`read_groups[0]`). If no processing is selected, each found @RG line will be its own entry in output array `read_groups`.",
            common: true
        }
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        Boolean format_for_star = true
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        if ~{format_for_star}; then
            samtools view -H ~{bam} \
                | grep "@RG" \
                | cut -f 2- \
                | sed -e 's/\t/ /g' \
                | awk '{print}' ORS=' , ' \
                | sed 's/ , $//' > read_groups.txt
        else
            samtools view -H ~{bam} | grep "@RG" > read_groups.txt
        fi
    >>>

    output {
        Array[String] read_groups = read_lines("read_groups.txt")
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'
        maxRetries: 1
    }
}

task split_string {
    # Currently (v1.1) no way to do this using the WDL standard library.
    # Revisit task in future version updates, can hopefully be replaced.
    meta {
        description: "Split a string into an array of strings based on a delimiter"
        outputs: {
            split_string: "Split string as an array"
        }
    }

    parameter_meta {
        input_string: "String to split on occurences of `delimiter`"
        delimiter: {
            description: "Delimiter on which to split `input_string`",
            common: true
        }
    }

    input {
        String input_string
        String delimiter = " , "
    }

    command <<<
        set -euo pipefail

        echo ~{input_string} | sed 's/~{delimiter}/\n/g' > split_strings.txt
    >>>

    output {
        Array[String] split_strings = read_lines("split_strings.txt")
    }

    runtime {
        memory: "4 GB"
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}

task calc_gene_lengths {
    meta {
        description: "Calculate gene lengths from a GTF feature file using the \"non-overlapping exonic length\" algorithm"
        help: "The \"non-overlapping exonic length\" algorithm can be implemented as the sum of each base covered by at least one exon; where each base is given a value of 1 regardless of how many exons overlap it."
        outputs: {
            gene_lengths: "A two column headered TSV file with gene names in the first column and feature lengths (as integers) in the second column"
        }
    }

    parameter_meta {
        gtf: "GTF feature file"
        outfile_name: "Name of the gene lengths file"
        idattr: {
            description: "GTF attribute to be used as feature ID. The value of this attribute will be used as the first column in the output file.",
            common: true
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
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/gtfparse:1.2.1--pyh864c0ab_0'
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
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'
        maxRetries: 1
    }
}

task add_to_bam_header {
    # TODO consider moving to samtools.wdl
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
        # TODO isn't there a builtin samtools command for this?
        samtools view -H ~{bam} > header.sam
        echo "~{additional_header}" >> header.sam
        samtools reheader -P header.sam ~{bam} > ~{outfile_name}
    >>>

    output {
        File reheadered_bam = outfile_name
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0'
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
        disk: "~{disk_size_gb} GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}

task make_coverage_regions_beds {
    # TODO make this customizable
    meta {
        description: "Takes in a GTF file, converts it to BED, then filters it down to two 3 column BED files: one of only 'exons', one of only 'CDS' regions"
        outputs: {
            bed: "Input GTF converted into BED format using the `gtf2bed` program",
            exon_bed: "3 column BED file corresponding to all 'exons' found in the input GTF",
            CDS_bed: "3 column BED file corresponding to all 'CDS' regions found in the input GTF"
        }
    }

    parameter_meta {
        gtf: "GTF feature file from which to derive coverage regions BED files"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File gtf
        Int modify_disk_size_gb = 0
    }

    Float gtf_size = size(gtf, "GiB")
    Int disk_size_gb = ceil(gtf_size * 2) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail

        BED=$(basename ~{gtf} '.gz').bed
        gunzip -c ~{gtf} | gtf2bed > "$BED"

        EXON=$(basename ~{gtf} '.gz').exon.bed
        awk '/\texon\t/ {print $1 "\t" $2 "\t" $3}' "$BED" > "$EXON"

        CDS=$(basename ~{gtf} '.gz').CDS.bed
        awk '/\tCDS\t/ {print $1 "\t" $2 "\t" $3}' "$BED" > "$CDS"
    >>>

    output {
        File bed = basename(gtf, '.gz') + ".bed"  # TODO I added this but I'm not sure I should have. It doesn't work with the task name
        File exon_bed = basename(gtf, '.gz') + ".exon.bed"
        File CDS_bed = basename(gtf, '.gz') + ".CDS.bed"
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/bedops:2.4.41--h9f5acd7_0'
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

    command <<<
        set -euo pipefail

        BAM="~{bam}" PREFIX="~{prefix}" FAST_MODE=~{fast_mode} python3 - <<END
import os  # lint-check: ignore
from collections import defaultdict  # lint-check: ignore

import pysam  # lint-check: ignore

bam_path = os.environ["BAM"]
bam = pysam.AlignmentFile(bam_path, "rb")

fast_mode = os.environ["FAST_MODE"] == "true"

if not fast_mode:
    tot_quals = defaultdict(lambda: 0)
    mapped_quals = defaultdict(lambda: 0)
    unmapped_quals = defaultdict(lambda: 0)
first_tot_quals = defaultdict(lambda: 0)
first_mapped_quals = defaultdict(lambda: 0)
first_unmapped_quals = defaultdict(lambda: 0)
middle_tot_quals = defaultdict(lambda: 0)
middle_mapped_quals = defaultdict(lambda: 0)
middle_unmapped_quals = defaultdict(lambda: 0)
last_tot_quals = defaultdict(lambda: 0)
last_mapped_quals = defaultdict(lambda: 0)
last_unmapped_quals = defaultdict(lambda: 0)
for read in bam:
    # only count primary alignments and unmapped reads
    if (read.is_secondary or read.is_supplementary) and not read.is_unmapped:
        continue

    cur_quals = read.query_qualities  # array of phred scores
    if not fast_mode:
        for qual in cur_quals:
            tot_quals[qual] += 1
            if read.is_unmapped:
                unmapped_quals[qual] += 1
            else:
                mapped_quals[qual] += 1

    first_score = cur_quals[0]
    first_tot_quals[first_score] += 1

    middle_pos = len(cur_quals) // 2  # middle base of read
    middle_score = cur_quals[middle_pos]
    middle_tot_quals[middle_score] += 1

    last_score = cur_quals[-1]
    last_tot_quals[last_score] += 1

    if read.is_unmapped:
        first_unmapped_quals[first_score] += 1
        middle_unmapped_quals[middle_score] += 1
        last_unmapped_quals[last_score] += 1
    else:
        first_mapped_quals[first_score] += 1
        middle_mapped_quals[middle_score] += 1
        last_mapped_quals[last_score] += 1

prefix = os.environ["PREFIX"]
outfile = open(prefix + ".global_PHRED_scores.tsv", "w")

# print header
header = ["sample"]
if not fast_mode:
    header += [
        "total average",
        "total median",
        "total stdev",
        "mapped average",
        "mapped median",
        "mapped stdev",
        "unmapped average",
        "unmapped median",
        "unmapped stdev",
    ]
header += [
    "first position total average",
    "first position total median",
    "first position total stdev",
    "first position mapped average",
    "first position mapped median",
    "first position mapped stdev",
    "first position unmapped average",
    "first position unmapped median",
    "first position unmapped stdev",
    "middle position total average",
    "middle position total median",
    "middle position total stdev",
    "middle position mapped average",
    "middle position mapped median",
    "middle position mapped stdev",
    "middle position unmapped average",
    "middle position unmapped median",
    "middle position unmapped stdev",
    "last position total average",
    "last position total median",
    "last position total stdev",
    "last position mapped average",
    "last position mapped median",
    "last position mapped stdev",
    "last position unmapped average",
    "last position unmapped median",
    "last position unmapped stdev",
]
print(
    "\t".join(header),
    file=outfile,
)
print(prefix, file=outfile, end="\t")


def stats_from_dict(score_dict):
    total_score = 0
    total_freq = 0
    freq_table = []
    for score, freq in score_dict.items():
        total_score += score * freq
        total_freq += freq
        freq_table.append((score, freq))

    if total_freq == 0:
        return -1, -1, -1
    avg = total_score / total_freq

    freq_table.sort(key=lambda entry: entry[0])

    median_pos = (total_freq + 1) / 2
    cumul_freq = 0
    median_found = False
    sum_freq_times_score_sqrd = 0
    for score, freq in freq_table:
        cumul_freq += freq
        if cumul_freq >= median_pos:
            if not median_found:
                median = score
            median_found = True
        sum_freq_times_score_sqrd += freq * (score**2)

    stdev = ((sum_freq_times_score_sqrd / total_freq) - (avg**2)) ** 0.5

    return avg, median, stdev


if not fast_mode:
    tot_avg, tot_median, tot_stdev = stats_from_dict(tot_quals)
    print(f"{tot_avg}", file=outfile, end="\t")
    print(f"{tot_median}", file=outfile, end="\t")
    print(f"{tot_stdev}", file=outfile, end="\t")

    mapped_avg, mapped_median, mapped_stdev = stats_from_dict(mapped_quals)
    print(f"{mapped_avg}", file=outfile, end="\t")
    print(f"{mapped_median}", file=outfile, end="\t")
    print(f"{mapped_stdev}", file=outfile, end="\t")

    unmapped_avg, unmapped_median, unmapped_stdev = stats_from_dict(unmapped_quals)
    print(f"{unmapped_avg}", file=outfile, end="\t")
    print(f"{unmapped_median}", file=outfile, end="\t")
    print(f"{unmapped_stdev}", file=outfile, end="\t")

first_tot_avg, first_tot_median, first_tot_stdev = stats_from_dict(first_tot_quals)
print(f"{first_tot_avg}", file=outfile, end="\t")
print(f"{first_tot_median}", file=outfile, end="\t")
print(f"{first_tot_stdev}", file=outfile, end="\t")

first_mapped_avg, first_mapped_median, first_mapped_stdev = stats_from_dict(
    first_mapped_quals
)
print(f"{first_mapped_avg}", file=outfile, end="\t")
print(f"{first_mapped_median}", file=outfile, end="\t")
print(f"{first_mapped_stdev}", file=outfile, end="\t")

first_unmapped_avg, first_unmapped_median, first_unmapped_stdev = stats_from_dict(
    first_unmapped_quals
)
print(f"{first_unmapped_avg}", file=outfile, end="\t")
print(f"{first_unmapped_median}", file=outfile, end="\t")
print(f"{first_unmapped_stdev}", file=outfile, end="\t")

middle_tot_avg, middle_tot_median, middle_tot_stdev = stats_from_dict(middle_tot_quals)
print(f"{middle_tot_avg}", file=outfile, end="\t")
print(f"{middle_tot_median}", file=outfile, end="\t")
print(f"{middle_tot_stdev}", file=outfile, end="\t")

middle_mapped_avg, middle_mapped_median, middle_mapped_stdev = stats_from_dict(
    middle_mapped_quals
)
print(f"{middle_mapped_avg}", file=outfile, end="\t")
print(f"{middle_mapped_median}", file=outfile, end="\t")
print(f"{middle_mapped_stdev}", file=outfile, end="\t")

middle_unmapped_avg, middle_unmapped_median, middle_unmapped_stdev = stats_from_dict(
    middle_unmapped_quals
)
print(f"{middle_unmapped_avg}", file=outfile, end="\t")
print(f"{middle_unmapped_median}", file=outfile, end="\t")
print(f"{middle_unmapped_stdev}", file=outfile, end="\t")

last_tot_avg, last_tot_median, last_tot_stdev = stats_from_dict(last_tot_quals)
print(f"{last_tot_avg}", file=outfile, end="\t")
print(f"{last_tot_median}", file=outfile, end="\t")
print(f"{last_tot_stdev}", file=outfile, end="\t")

last_mapped_avg, last_mapped_median, last_mapped_stdev = stats_from_dict(
    last_mapped_quals
)
print(f"{last_mapped_avg}", file=outfile, end="\t")
print(f"{last_mapped_median}", file=outfile, end="\t")
print(f"{last_mapped_stdev}", file=outfile, end="\t")

last_unmapped_avg, last_unmapped_median, last_unmapped_stdev = stats_from_dict(
    last_unmapped_quals
)
print(f"{last_unmapped_avg}", file=outfile, end="\t")
print(f"{last_unmapped_median}", file=outfile, end="\t")
print(f"{last_unmapped_stdev}", file=outfile)  # end="\n"

outfile.close()

END
    >>>

    output {
        File phred_scores = "~{outfile_name}"
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}

task qc_summary {
    # TODO not sure why I implemented as a JSON. Wouldn't TSV be easier to work with? Talk to Delaram/David
    #   Delaram+David okayed a switch to TSV
    meta {
        description: "**[OUT OF DATE]** This WDL task pulls out keys metrics that can provide a high level overview of the sample, without needing to examine the entire MultiQC report. Currently, these key metrics come from Qualimap and ngsderive."
    }

    input {
        File multiqc_tar_gz
        String outfile_name = basename(multiqc_tar_gz, ".multiqc.tar.gz") + ".qc_summary.json"
    }

    String sample_name = basename(multiqc_tar_gz, ".multiqc.tar.gz")

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
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}
