import pysam
import json


def main(bam_path, outfile_path):
    sam = pysam.AlignmentFile(bam_path, "rb")

    out_file = open(outfile_path, "w")
    header = sam.header.to_dict()["RG"]
    modified_header = []
    for read_group in sorted(header, key=lambda d: d["ID"]):
        modified_header.append(
            {k: v.upper() if k == "PL" else v for k, v in read_group.items()}
        )
    json.dump(modified_header, out_file)
    out_file.close()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bam_path",
        type=str,
    )
    parser.add_argument(
        "outfile_path",
        type=str,
    )
    args = parser.parse_args()

    main(args.bam_path, args.outfile_path)
