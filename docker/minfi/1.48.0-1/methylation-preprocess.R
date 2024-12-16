options(error = function() traceback(3))

library(minfi)
library(dplyr)
library(data.table)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify options
parser$add_argument("--idat_base", type="character", help="Base path to the idat files")
parser$add_argument("--out_base", type="character", help="Base path to the output files")
# parse the arguments
args <- parser$parse_args()

# These arrays have two types of probes (helpfully named Infinium type 1 and Infinium type 2).
# The probes also cover a varying number of CpG sites.
# The normalization method works by taking N probes of each type that each have 1, 2, or 3 CpGs.
# So it selects 6N probes in total. This selection is "random" if the seed is not fixed, of course.
# Since we're processing samples individually in this step instead of as a cohort,
# we want the seed to be consistent so the same set of probes is chosen for each sample.
set.seed(1)

dir=getwd()

RGSet <- read.metharray(args$idat_base, verbose = TRUE, force = TRUE)
saveRDS(RGSet, paste0(args$out_base, ".RGSet.rds"))

# The manifest is needed by preprocessRAW
manifest <- getManifest(RGSet)
manifest

# Load raw data into a MethylSet object be converting red/green
# channels into a matrix of methlyated and unmethylated signals.
MSet <- preprocessRaw(RGSet)
saveRDS(MSet, paste0(args$out_base, ".MSet.rds"))

# Convert to a RatioSet
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
saveRDS(RSet, paste0(args$out_base, ".RSet.rds"))

# Add genomic coordinates to each probe (plus additional annotation)
GRset <- mapToGenome(RSet)
saveRDS(GRset, paste0(args$out_base, ".GRSet.rds"))

# Take the genomic mapped RatioSet and fill Beta values (non-normalized).
# Get the NON-normalized beta values:
beta <- getBeta(GRset)
write.csv(beta, paste0(args$out_base, ".beta.csv"))

# Get M-value matrix and copy-number matrix
# Get M and CN vals if interested:
M <- getM(GRset)
write.csv(M, paste0(args$out_base, ".m_values.csv"))
CN <- getCN(GRset)
write.csv(CN, paste0(args$out_base, ".cn_values.csv"))

# Get sample names and probe names
sampleNames <- sampleNames(GRset)
write.csv(sampleNames, paste0(args$out_base, ".sampleNames.csv"))
probeNames <- featureNames(GRset)
write.csv(probeNames, paste0(args$out_base, ".probeNames.csv"))

gr <- granges(GRset)
write.csv(gr, paste0(args$out_base, ".gr.csv"))

annotation <- getAnnotation(GRset)
write.csv(annotation, paste0(args$out_base, ".annotation.csv"))

# Perform SWAN normalization on beta values
GRset.swan_norm <- preprocessSWAN(RGSet)
write.csv(GRset.swan_norm, paste0(args$out_base, ".GRset.swan_norm.csv"))
beta_swan_norm <- getBeta(GRset.swan_norm)

# Write the normalized beta-values that have NOT yet had
# low-variance probes filtered out
write.csv(beta_swan_norm, paste0(args$out_base, ".beta_swan_norm_unfiltered.csv"))

# Write the normalized beta-values that have NOT yet had
# low-variance probes filtered out
# Filter to only those that are mappable to the genome.
RSet <- ratioConvert(GRset.swan_norm)
GRset <- mapToGenome(RSet)
beta_swan_norm <- getBeta(GRset)
colnames(beta_swan_norm) <- colnames(CN)
write.csv(beta_swan_norm[order(rownames(beta_swan_norm)),], paste0(args$out_base, ".beta_swan_norm_unfiltered.genomic.csv"))