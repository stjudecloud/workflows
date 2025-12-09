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

rg_set <- read.metharray(args$idat_base, verbose = TRUE, force = TRUE)
saveRDS(rg_set, paste0(args$out_base, ".RGSet.rds"))

# Generate detection p-values for all probed positions
det_p <- detectionP(rg_set)
det_p <-
  det_p[order(rownames(det_p)), , drop = FALSE]

# The manifest is needed by preprocessRAW
manifest <- getManifest(rg_set)
manifest

# Load raw data into a MethylSet object be converting red/green
# channels into a matrix of methlyated and unmethylated signals.
m_set <- preprocessRaw(rg_set)
saveRDS(m_set, paste0(args$out_base, ".MSet.rds"))

# Convert to a RatioSet
r_set <- ratioConvert(m_set, what = "both", keepCN = TRUE)
saveRDS(r_set, paste0(args$out_base, ".RSet.rds"))

# Add genomic coordinates to each probe (plus additional annotation)
gr_set <- mapToGenome(r_set)
saveRDS(gr_set, paste0(args$out_base, ".GRSet.rds"))

# Get the list of sites with SNPs
snps <- getSnpInfo(gr_set)
gr_set_snps <- addSnpInfo(gr_set)
# Remove probes with SNPs at CpG or single-base extension sites
gr_set_snps <- dropLociWithSnps(gr_set_snps, snps=c("SBE","CpG"), maf=0)
probes_without_snps <- featureNames(gr_set_snps)

# Take the genomic mapped RatioSet and fill Beta values (non-normalized).
# Get the NON-normalized beta values:
beta <- getBeta(gr_set)
write.csv(beta, paste0(args$out_base, ".beta.csv"))

# Get M-value matrix and copy-number matrix
# Get M and CN vals if interested:
m <- getM(gr_set)
write.csv(m, paste0(args$out_base, ".m_values.csv"))
cn <- getCN(gr_set)
write.csv(cn, paste0(args$out_base, ".cn_values.csv"))

# Get sample names and probe names
sample_names <- sampleNames(gr_set)
write.csv(sample_names, paste0(args$out_base, ".sampleNames.csv"))
probe_names <- featureNames(gr_set)
write.csv(probe_names, paste0(args$out_base, ".probeNames.csv"))

# Write probes with SNPs
probes_with_snps <- setdiff(probe_names, probes_without_snps)
write.csv(probes_with_snps, paste0(args$out_base, ".probes_with_snps.csv"))
write.csv(
  probes_without_snps,
  paste0(args$out_base, ".probes_without_snps.csv")
)

gr <- granges(gr_set)
write.csv(gr, paste0(args$out_base, ".gr.csv"))

annotation <- getAnnotation(gr_set)
write.csv(annotation, paste0(args$out_base, ".annotation.csv"))

# Perform SWAN normalization on beta values
gr_set_swan_norm <- preprocessSWAN(rg_set)
beta_swan_norm <- getBeta(gr_set_swan_norm)

# Write the normalized beta-values that have NOT yet had
# low-variance probes filtered out
write.csv(
  beta_swan_norm,
  paste0(args$out_base, ".beta_swan_norm_unfiltered.csv")
)

# Write the normalized beta-values that have NOT yet had
# low-variance probes filtered out
# Filter to only those that are mappable to the genome.
r_set <- ratioConvert(gr_set_swan_norm)
gr_set <- mapToGenome(r_set)
beta_swan_norm <- getBeta(gr_set)
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
