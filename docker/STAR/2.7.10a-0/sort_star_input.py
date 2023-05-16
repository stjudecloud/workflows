"""Sort FastQ pairs and read group information prior to being passed to STAR.

We discovered that when providing STAR with read group split FastQs and also
supplying STAR with header information via the `--outSAMattrRGline` parameter,
that the read group information for individual reads may become jumbled.

The simplest solution for this was to first pass exactly what would be passed
to STAR into this program, sort it correctly, and write the reordered
arguments out to files which can be supplied to STAR via the `cat` command.
"""

import argparse
import os
from typing import List, Tuple, Any


def sort_fastqs(fastq_files: List[str]) -> List[str]:
    """Sort and return a list of file paths according to their basename."""
    fastq_basenames = [fastq.split(os.sep)[-1] for fastq in fastq_files]

    sorted_filepaths, sorted_basenames = zip(
        *sorted(zip(fastq_files, fastq_basenames), key=lambda item: item[1])
    )

    return sorted_filepaths


def sort_read_groups(read_groups_string: str) -> Tuple[List[str], List[str]]:
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
    rgids = [flags.split(" ")[0].split(":")[1] for flags in read_groups]

    sorted_read_groups, sorted_rgids = zip(
        *sorted(zip(read_groups, rgids), key=lambda item: item[1])
    )

    return (sorted_read_groups, sorted_rgids)


def validate(
    read_one_fastqs: List[str], read_two_fastqs: List[str], rgids: List[str]
) -> None:
    """Ensure that the final strings are ready to be passed to STAR.

    The first check is that there is the same number of read groups as
    FastQ pairs. Next, it's checked that each read group ID is
    present in a synced FastQ pair.

    Args:
        read_one_fastqs (list): list of file paths
        read_two_fastqs (list): list of file paths
        rgids (list): list of values to the 'ID' key
    """
    if len(read_one_fastqs) != len(rgids):
        raise argparse.ArgumentError(
            "Must have same number of read groups as FastQ pairs"
        )

    for i, id in enumerate(rgids):
        if (id not in read_one_fastqs[i]) or (id not in read_two_fastqs[i]):
            raise SystemExit(
                "Error: There's a mismatch between "
                "read group IDs and FastQ file names"
            )


def write_outfiles(
    read_one_fastqs: List[str], read_two_fastqs: List[str], read_groups: List[str]
) -> None:
    """Concatenate each element in each list and write results to files."""
    read_one_file = open("read_one_fastqs_sorted.txt", "w")
    read_one_file.write(",".join(read_one_fastqs))
    read_one_file.close()

    read_two_file = open("read_two_fastqs_sorted.txt", "w")
    read_two_file.write(",".join(read_two_fastqs))
    read_two_file.close()

    read_group_file = open("read_groups_sorted.txt", "w")
    read_group_file.write(" , ".join(read_groups))
    read_group_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--read_one_fastqs",
        required=True,
        help="Comma delimited (without spaces) list of read one FastQ file paths",
    )
    parser.add_argument(
        "--read_two_fastqs",
        required=False,
        help="Comma delimited (without spaces) list of read two FastQ file paths",
    )
    parser.add_argument(
        "--read_groups",
        required=False,
        help="Exactly what would be passed to STAR "
        "after `--outSAMattrRGline`.\nSince this will have spaces, "
        "and this is a single argument, quote the entire string.",
    )

    args = parser.parse_args()

    read_one_fastqs = args.read_one_fastqs.split(",")
    if args.read_two_fastqs:
        read_two_fastqs = args.read_two_fastqs.split(",")

    if args.read_two_fastqs:
        if len(read_one_fastqs) != len(read_two_fastqs):
            raise argparse.ArgumentError(
                "Must have the same number of read one FastQs as read two FastQs"
            )

    sorted_read_one_fastqs = sort_fastqs(read_one_fastqs)
    sorted_read_two_fastqs = [""]
    if args.read_two_fastqs:
        sorted_read_two_fastqs = sort_fastqs(read_two_fastqs)

    sorted_read_groups = ["-"]  # Indicates default to STAR
    if args.read_groups:
        sorted_read_groups, sorted_rgids = sort_read_groups(args.read_groups)
        validate(sorted_read_one_fastqs, sorted_read_two_fastqs, sorted_rgids)

    write_outfiles(sorted_read_one_fastqs, sorted_read_two_fastqs, sorted_read_groups)
