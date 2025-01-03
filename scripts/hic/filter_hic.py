
from collections import defaultdict
import argparse

def get_args():
    parser = argparse.ArgumentParser(
        description="Filter Hi-C data.")
    parser.add_argument(
        "--prefix", type=str, help="Prefix for output file.")
    parser.add_argument(
        "--all_valid_pairs", type=str, help="All valid pairs file.")
    parser.add_argument(
        "--filter_pairs", type=str, help="Filter pairs file.")
    
    args = parser.parse_args()
    return args
    
if __name__ == "__main__":
    args = get_args()
    
    f=open(args.filter_pairs)
    blackID=defaultdict(int)
    while True:
        line=f.readline()
        if not line:
            break
        cols=line.strip().split("\t")
        ID=cols[0]
        blackID[ID]=1

    f2=open(args.all_valid_pairs)
    outfile1=args.prefix + ".allValidPairs.filtered"
    outfile2=args.prefix + ".allValidPairs.removed"
    of1=open(outfile1,'w')
    of2=open(outfile2,'w')
    while True:
        line=f2.readline()
        if not line:
            break
        cols=line.strip().split("\t")
        id=cols[0]
        if id in blackID:
            of2.write(line)
        else:
            of1.write(line)