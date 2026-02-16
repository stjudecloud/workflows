library(minfi)
library(dplyr)
library(data.table)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

suppressPackageStartupMessages(library("argparse"))
# create parser object
parser <- ArgumentParser()
# specify options
parser$add_argument(
  "--idat_base",
  type = "character",
  help = "Base path to the idat files"
)
parser$add_argument(
  "--out_base",
  type = "character",
  help = "Base path to the output files"
)
parser$add_argument(
  "--seed",
  type = "numeric",
  default = 1,
  help = "Seed for random number generation"
)
# parse the arguments
args <- parser$parse_args()

# These arrays have two types of probes (named Infinium type 1 and Infinium
# type 2). The probes also cover a varying number of CpG sites. The
# normalization method works by taking N probes of each type that each have 1,
# 2, or 3 CpGs. So it selects 6N probes in total. This selection is "random"
# if the seed is not fixed, of course. Since we're processing samples
# individually in this step instead of as a cohort. We want the seed to be
# consistent so the same set of probes is chosen for each sample.
set.seed(args$seed)

dir <- getwd()

red_green_channel_set <-
  read.metharray(args$idat_base, verbose = TRUE, force = TRUE)
saveRDS(red_green_channel_set, paste0(args$out_base, ".RGSet.rds"))

# Generate detection p-values for all probed positions
det_p <- detectionP(red_green_channel_set)
det_p <-
  det_p[order(rownames(det_p)), , drop = FALSE]

# The manifest is needed by preprocessRAW
manifest <- getManifest(red_green_channel_set)
manifest

# Load raw data into a MethylSet object be converting red/green
# channels into a matrix of methlyated and unmethylated signals.
methyl_set <- preprocessRaw(red_green_channel_set)
saveRDS(methyl_set, paste0(args$out_base, ".MSet.rds"))

# Convert to a RatioSet
ratio_set <- ratioConvert(methyl_set, what = "both", keepCN = TRUE)
saveRDS(ratio_set, paste0(args$out_base, ".RSet.rds"))

# Add genomic coordinates to each probe (plus additional annotation)
genomic_ratio_set <- mapToGenome(ratio_set)
saveRDS(genomic_ratio_set, paste0(args$out_base, ".GRSet.rds"))

non_genomic_probes <-
  setdiff(featureNames(ratio_set), featureNames(genomic_ratio_set))

# Get the list of sites with SNPs
snps <- getSnpInfo(genomic_ratio_set)
genomic_ratio_set_snps <- addSnpInfo(genomic_ratio_set)
# Remove probes with SNPs at CpG or single-base extension sites
genomic_ratio_set_snps <-
  dropLociWithSnps(genomic_ratio_set_snps, snps = c("SBE", "CpG"), maf = 0)
probes_without_snps <- featureNames(genomic_ratio_set_snps)

# Take the genomic mapped RatioSet and fill Beta values (non-normalized).
# Get the NON-normalized beta values:
beta <- getBeta(genomic_ratio_set)
write.csv(beta, paste0(args$out_base, ".beta.csv"))

# Get M-value matrix and copy-number matrix
# Get M and CN vals if interested:
m <- getM(genomic_ratio_set)
write.csv(m, paste0(args$out_base, ".m_values.csv"))
cn <- getCN(genomic_ratio_set)
write.csv(cn, paste0(args$out_base, ".cn_values.csv"))

# Get sample names and probe names
sample_names <- sampleNames(genomic_ratio_set)
write.csv(sample_names, paste0(args$out_base, ".sampleNames.csv"))
probe_names <- featureNames(genomic_ratio_set)
write.csv(probe_names, paste0(args$out_base, ".probeNames.csv"))

# Write probes with SNPs
probes_with_snps <- setdiff(probe_names, probes_without_snps)
write.table(
  probes_with_snps,
  paste0(args$out_base, ".probes_with_snps.tab"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
)
write.table(
  probes_without_snps,
  paste0(args$out_base, ".probes_without_snps.tab"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
)

# Write non-genomic probe list
write.table(
  non_genomic_probes,
  paste0(args$out_base, ".non_genomic_probes.tab"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
)

gr <- granges(genomic_ratio_set)
write.csv(gr, paste0(args$out_base, ".gr.csv"))

annotation <- getAnnotation(genomic_ratio_set)
write.csv(annotation, paste0(args$out_base, ".annotation.csv"))

# Perform SWAN normalization on beta values
genomic_ratio_set_swan_norm <- preprocessSWAN(red_green_channel_set)
beta_swan_norm <- getBeta(genomic_ratio_set_swan_norm)

# Write the normalized beta-values that have NOT yet had
# low-variance probes filtered out
write.csv(
  beta_swan_norm,
  paste0(args$out_base, ".beta_swan_norm_unfiltered.csv")
)

# Write the normalized beta-values that have NOT yet had
# low-variance probes filtered out
# Filter to only those that are mappable to the genome.
ratio_set <- ratioConvert(genomic_ratio_set_swan_norm)
genomic_ratio_set <- mapToGenome(ratio_set)
beta_swan_norm <- getBeta(genomic_ratio_set)
beta_swan_norm <-
  beta_swan_norm[order(rownames(beta_swan_norm)), , drop = FALSE]
write.csv(
  beta_swan_norm,
  paste0(args$out_base, ".beta_swan_norm_unfiltered.genomic.csv")
)

genomic_probes <- rownames(beta_swan_norm)
all_probes <- rownames(det_p)
non_genomic_probes <- setdiff(all_probes, genomic_probes)

# Filter probe p-values to only those with genomic coordinates
det_p <- det_p[!(row.names(det_p) %in% non_genomic_probes), , drop = FALSE]
det_p <-
  det_p[order(rownames(det_p)), , drop = FALSE]

write.csv(det_p, paste0(args$out_base, ".detectionP.csv"))
