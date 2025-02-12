from collections import defaultdict

import pysam


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
    median = None
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


def main(bam_path, prefix, fast_mode):
    bam = pysam.AlignmentFile(bam_path, "rb")

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

    middle_tot_avg, middle_tot_median, middle_tot_stdev = stats_from_dict(
        middle_tot_quals
    )
    print(f"{middle_tot_avg}", file=outfile, end="\t")
    print(f"{middle_tot_median}", file=outfile, end="\t")
    print(f"{middle_tot_stdev}", file=outfile, end="\t")

    middle_mapped_avg, middle_mapped_median, middle_mapped_stdev = stats_from_dict(
        middle_mapped_quals
    )
    print(f"{middle_mapped_avg}", file=outfile, end="\t")
    print(f"{middle_mapped_median}", file=outfile, end="\t")
    print(f"{middle_mapped_stdev}", file=outfile, end="\t")

    (
        middle_unmapped_avg,
        middle_unmapped_median,
        middle_unmapped_stdev,
    ) = stats_from_dict(middle_unmapped_quals)
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


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bam_path",
        type=str,
        help="Path to the BAM file to calculate global PHRED scores for",
    )
    parser.add_argument(
        "prefix",
        type=str,
        help="Prefix for the output file",
    )
    parser.add_argument(
        "--fast_mode",
        action="store_true",
        help="Only calculate the average, median, and standard deviation of the global PHRED scores",
    )
    args = parser.parse_args()

    main(args.bam_path, args.prefix, args.fast_mode)
