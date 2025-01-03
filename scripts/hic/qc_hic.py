import os
import argparse

def get_args():
    parser = argparse.ArgumentParser(
        description="QC a Hi-C experiment.")
    parser.add_argument(
        "--all_valid_pairs_stats", type=str, help="All valid pairs stats file.")
    parser.add_argument(
        "--pairing_stats", type=str, help="Pairing stats file.")
    parser.add_argument(
        "--mapping_stats_read1", type=str, help="Mapping stats for read 1.")
    parser.add_argument(
        "--mapping_stats_read2", type=str, help="Mapping stats for read 2.")
    parser.add_argument(
        "--fithichip_bed", type=str, help="FithiChIP bed file.")
    parser.add_argument(
        "--fithichip_q01_bed", type=str, help="FithiChIP q01 bed file.")
    parser.add_argument(
        "--peaks_bed", type=str, help="Peaks bed file.")
    parser.add_argument(
        "--prefix", type=str, help="Prefix for output file.")
        
    
    args = parser.parse_args()
    
    return args

if __name__ == "__main__":
    args = get_args()
    FILES = [args.all_valid_pairs_stats, args.pairing_stats, args.mapping_stats_read1, args.mapping_stats_read2]
    VARIABLE = ["valid_interaction", "cis_shortRange", "cis_longRange", "Total_pairs_processed", "Reported_pairs", "total_R2", "mapped_R2", "total_R1", "mapped_R1"]
    RESULTS = {}
    for eachfile in FILES:
        with open(eachfile) as f:
            for line in f:
                array = line.split()
                if array[0] in VARIABLE:
                    RESULTS[array[0]] = int(array[1])

    PERCENTAGES = {}
    COUNTS = {}
    VERDICT = {}
    REP_VAR = ["R1_aligned", "R2_aligned", "valid_interactionPairs", "cis_shortRange", "cis_longRange"]

    if os.path.isfile(args.fithichip_bed):
        with open(args.fithichip_q01_bed) as fithichip:
            LOOPS_SIGNIFICANT = len(fithichip.readlines())-1
        with open(args.fithichip_bed) as fithichip:
            LOOPS = len(fithichip.readlines())-1
    if os.path.isfile(args.peaks_bed):
        with open(args.peaks_bed) as peaks:
            PEAKS = len(peaks.readlines())-1

    PERCENTAGES["R1_aligned"] = round((RESULTS["mapped_R1"]*100)/RESULTS["total_R1"])
    PERCENTAGES["R2_aligned"] = round((RESULTS["mapped_R2"]*100)/RESULTS["total_R2"])
    PERCENTAGES["valid_interactionPairs"] = round((RESULTS["valid_interaction"]*100)/RESULTS["Reported_pairs"])
    PERCENTAGES["cis_shortRange"] = round((RESULTS["cis_shortRange"]*100)/RESULTS["valid_interaction"])
    PERCENTAGES["cis_longRange"] = round((RESULTS["cis_longRange"]*100)/RESULTS["valid_interaction"])
    COUNTS["R1_aligned"] = RESULTS["mapped_R1"]
    COUNTS["R2_aligned"] = RESULTS["mapped_R2"]
    COUNTS["valid_interactionPairs"] = RESULTS["valid_interaction"]
    COUNTS["cis_shortRange"] = RESULTS["cis_shortRange"]
    COUNTS["cis_longRange"] = RESULTS["cis_longRange"]

    for aligned in ["R1_aligned", "R2_aligned"]:
        if PERCENTAGES[aligned] > 80:
            VERDICT[aligned] = "GOOD"
        else:
            VERDICT[aligned] = "BAD"

    if PERCENTAGES["valid_interactionPairs"] > 50:
        VERDICT["valid_interactionPairs"] = "GOOD"
    else:
        VERDICT["valid_interactionPairs"] = "BAD"

    if PERCENTAGES["cis_shortRange"] > 50:
        VERDICT["cis_shortRange"] = "BAD"
    elif PERCENTAGES["cis_shortRange"] > 30:
        VERDICT["cis_shortRange"] = "MARGINAL"
    else:
        VERDICT["cis_shortRange"] = "GOOD"

    if PERCENTAGES["cis_longRange"] > 40:
        VERDICT["cis_longRange"] = "GOOD"
    elif PERCENTAGES["cis_longRange"] > 20:
        VERDICT["cis_longRange"] = "MARGINAL"
    else:
        VERDICT["cis_longRange"] = "BAD"

    REPORT = open(args.prefix + "_QCreport.txt", 'w')
    REPORT.write("STAT\tCOUNTS\tPERCENTAGE\tVERDICT\n")
    REPORT.write("Total_pairs_processed\t" + str(RESULTS["Total_pairs_processed"]) + "\n")
    for var in REP_VAR:
        REPORT.write(var + "\t" + str(COUNTS[var]) + "\t" + str(PERCENTAGES[var]) + "\t" + VERDICT[var] + '\n')

    if os.path.isfile(args.peaks_bed):
        REPORT.write("peaks\t" + str(PEAKS) + "\n")
    if os.path.isfile(args.fithichip_bed):
        REPORT.write("loops\t" + str(LOOPS) + "\n")
        REPORT.write("loops_significant\t" + str(LOOPS_SIGNIFICANT) + "\n")