"""Check FASTQ pairs and read group information prior to further processing.

This script checks various assumptions we make about FASTQ and RG record concordance.
- each read1 FASTQ must correspond to exactly one RG record
    - this correspondence is encoded two ways:
       1. they share the same index of their respective arrays
       2. the `ID` field of the RG record is present somewhere in the read1 FASTQ basename
- if read2 FASTQs are present each must correspond to exactly one read1 FASTQ/RG record pair
    - this correspondence is encoded the same way as it is for read1 FASTQs
- the `ID` field of each RG record must be the first field
- every `ID` must be unique

Read group records may be optionally prefixed with `@RG` and may use either tabs or spaces as delimiters.
"""

import argparse
import os
from typing import List


def validate(
    read1_fastqs: List[str],
    read2_fastqs: List[str],
    rgids: List[str],
) -> None:
    """Ensure that each FASTQ basename contains the expected RG ID and that each ID is unique.

    Assumes that read1_fastqs and rgids are the same length
    and that read2_fastqs is empty or the same length as the other lists.
    """
    checking_r2 = len(read2_fastqs) > 0
    for i, rgid in enumerate(rgids):
        if (rgid not in read1_fastqs[i]) or (
            checking_r2 and rgid not in read2_fastqs[i]
        ):
            raise SystemExit(f"{rgid} not found in corresponding FASTQ basename")
    
    if len(rgids) != len(set(rgids)):
        raise SystemExit("Each ID field must be unique")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--read-one-fastqs",
        required=True,
        help="Comma delimited list of read one FASTQ file paths. Entire argument should be quoted in the call to this script.",
    )
    parser.add_argument(
        "--read-two-fastqs",
        required=False,
        help="Comma delimited list of read two FASTQ file paths. Entire argument should be quoted in the call to this script.",
    )
    parser.add_argument(
        "--read-groups",
        required=True,
        help="Comma delimited list of read group records. Individual fields of each RG are whitespace delimited. Entire argument should be quoted in the call to this script.",
    )

    args = parser.parse_args()

    read1_fastqs = [
        fq.strip().split(os.sep)[-1] for fq in args.read_one_fastqs.split(",")
    ]
    read2_fastqs = []
    if args.read_two_fastqs:
        read2_fastqs = [
            fq.strip().split(os.sep)[-1] for fq in args.read_two_fastqs.split(",")
        ]

    if read2_fastqs != []:
        if len(read1_fastqs) != len(read2_fastqs):
            raise SystemExit(
                "Must have the same number of read one FASTQs as read two FASTQs"
            )

    rg_records = [rg.strip() for rg in args.read_groups.split(",")]

    if len(rg_records) != len(read1_fastqs):
        raise SystemExit("Must have the same number of RG records as read one FASTQs")

    rgids = []
    for rg in rg_records:
        split = rg.split()
        if split[0] == "@RG":
            first = split[1]
        else:
            first = split[0]
        field, id = first.split(":")
        if field != "ID":
            raise SystemExit("ID field must be the first field for each RG record")
        rgids.append(id)

    validate(read1_fastqs, read2_fastqs, rgids)
