"""Sort FASTQ pairs and read group information prior to being passed to STAR.

We discovered that when providing STAR with read group split FASTQs and also
supplying STAR with header information via the `--outSAMattrRGline` parameter,
that the read group information for individual reads may become jumbled.

The simplest solution for this was to first pass exactly what would be passed
to STAR into this program, sort it correctly, and write the reordered
arguments out to files which can be supplied to STAR via the `cat` command.
"""

import argparse
import os
from typing import List, Tuple


def sort_fastqs(fastq_files: List[str]) -> Tuple[str, ...]:
    """Sort and return a list of file paths according to their basename."""
    fastq_basenames = [fastq.split(os.sep)[-1] for fastq in fastq_files]

    sorted_filepaths, _sorted_basenames = zip(
        *sorted(zip(fastq_files, fastq_basenames), key=lambda item: item[1])
    )

    return sorted_filepaths


def sort_read_groups(
    read_groups_string: str,
) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    """Reorder the argument to the STAR `--outSAMattrRGline` parameter.

    The sort line could be substituted by using the `sorted()`
    method on both lists, since what's being sorted on comes
    first in the strings.
    However the `rgids` list is needed later in the program,
    and using the same method as `sort_fastqs()` contributes to consistency.

    Args:
        read_groups_string (str): Exactly what would be passed to STAR
            after `--outSAMattrRGline`.

    Returns:
        list: A list consisting of the input string split on ' , ',
            and then sorted by the read group ID values.
        list: A sorted list of values to the 'ID' key.
    """
    read_groups = [rg for rg in read_groups_string.split(" , ")]
    rgids = [
        flags.split(" ")[0].split(":")[1] for flags in read_groups
    ]  # TODO this breaks if flags are out of expected order

    sorted_read_groups, sorted_rgids = zip(
        *sorted(zip(read_groups, rgids), key=lambda item: item[1])
    )

    return (sorted_read_groups, sorted_rgids)


def validate(
    read_one_fastqs: Tuple[str, ...],
    read_two_fastqs: Tuple[str, ...],
    rgids: Tuple[str, ...],
) -> None:
    """Ensure that the final strings are ready to be passed to STAR.

    The first check is that there is the same number of read groups as
    FASTQ pairs. Next, it's checked that each read group ID is
    present in a synced FASTQ pair.

    Args:
        read_one_fastqs (list): list of file paths
        read_two_fastqs (list): list of file paths
        rgids (list): list of values to the 'ID' key
    """
    if len(read_one_fastqs) != len(rgids):
        raise SystemExit("Must have same number of read groups as FASTQ pairs")

    for i, rgid in enumerate(rgids):
        if (rgid not in read_one_fastqs[i]) or (rgid not in read_two_fastqs[i]):
            raise SystemExit(
                "Error: There's a mismatch between "
                + "read group IDs and FASTQ file names"
            )


def write_outfiles(
    read_one_fastqs: Tuple[str, ...],
    read_two_fastqs: Tuple[str, ...],
    read_groups: Tuple[str, ...],
) -> None:
    """Concatenate each element in each list and write results to files."""
    read_one_file = open("read_one_fastqs_sorted.txt", "w", encoding="utf-8")
    read_one_file.write(",".join(read_one_fastqs))
    read_one_file.write("\n")
    read_one_file.close()

    read_two_file = open("read_two_fastqs_sorted.txt", "w", encoding="utf-8")
    read_two_file.write(",".join(read_two_fastqs))
    read_two_file.write("\n")
    read_two_file.close()

    read_group_file = open("read_groups_sorted.txt", "w", encoding="utf-8")
    read_group_file.write(" , ".join(read_groups))
    read_group_file.write("\n")
    read_group_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--read-one-fastqs",
        required=True,
        help="Comma delimited (without spaces) list of read one FASTQ file paths",
    )
    parser.add_argument(
        "--read-two-fastqs",
        required=False,
        help="Comma delimited (without spaces) list of read two FASTQ file paths",
    )
    parser.add_argument(
        "--read-groups",
        required=False,
        help="Exactly what would be passed to STAR "
        "after `--outSAMattrRGline`.\nSince this will have spaces, "
        "and this is a single argument, quote the entire string.",
    )

    args = parser.parse_args()

    read1_fastqs = args.read_one_fastqs.split(",")
    if args.read_two_fastqs:
        read2_fastqs = args.read_two_fastqs.split(",")

    if args.read_two_fastqs:
        if len(read1_fastqs) != len(read2_fastqs):
            raise SystemExit(
                "Must have the same number of read one FASTQs as read two FASTQs"
            )

    sorted_read1_fastqs = sort_fastqs(read1_fastqs)
    sorted_read2_fastqs = tuple("")
    if args.read_two_fastqs:
        sorted_read2_fastqs = sort_fastqs(read2_fastqs)

    sorted_read_group_strs = tuple("-")  # Indicates default to STAR
    if args.read_groups:
        sorted_read_group_strs, sorted_read_group_ids = sort_read_groups(
            args.read_groups
        )
        validate(sorted_read1_fastqs, sorted_read2_fastqs, sorted_read_group_ids)

    write_outfiles(sorted_read1_fastqs, sorted_read2_fastqs, sorted_read_group_strs)
