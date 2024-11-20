version 1.1

workflow methylation_preprocess {
    meta {
        description: "Process raw IDAT files from the Illumina methylation array for a single sample"
        outputs: {
            beta_swan_norm_unfiltered: "Normalized beta values for all probes",
            beta_swan_norm_unfiltered_genomic: "Normalized beta values for all probes that map to the genome",
            mset: "MethylSet object",
            rgset: "RGSet object",
            rset: "RatioSet object",
            annotation: "Annotation object",
            beta: "Beta values",
            cn_values: "Copy number values",
            m_values: "M values",
            pheno: "Phenotype data",
            pheno_data: "Phenotype data",
            probe_names: "Probe names",
            sample_names: "Sample names",
        }
    }

    parameter_meta {
        idats: "Array of raw IDAT files from the Illumina methlyation array"
    }

    input {
        Pair[File, File] idats
    }

    call process_raw_idats { input:
        idats,
    }

    #@ except: LineWidth
    output {
        File beta_swan_norm_unfiltered = process_raw_idats.beta_swan_norm_unfiltered
        File beta_swan_norm_unfiltered_genomic = process_raw_idats.beta_swan_norm_unfiltered_genomic
        File mset = process_raw_idats.mset
        File rgset = process_raw_idats.rgset
        File rset = process_raw_idats.rset
        File annotation = process_raw_idats.annotation
        File beta = process_raw_idats.beta
        File cn_values = process_raw_idats.cn_values
        File m_values = process_raw_idats.m_values
        File pheno = process_raw_idats.pheno
        File pheno_data = process_raw_idats.pheno_data
        File probe_names = process_raw_idats.probe_names
        File sample_names = process_raw_idats.sample_names
    }
}

task process_raw_idats {
    meta {
        description: "Process raw IDAT files from the Illumina methylation array for a single sample"
        outputs: {
            beta_swan_norm_unfiltered: "Normalized beta values for all probes",
            beta_swan_norm_unfiltered_genomic: "Normalized beta values for all probes that map to the genome",
            mset: "MethylSet object",
            rgset: "RGSet object",
            rset: "RatioSet object",
            annotation: "Annotation object",
            beta: "Beta values",
            cn_values: "Copy number values",
            m_values: "M values",
            pheno: "Phenotype data",
            pheno_data: "Phenotype data",
            probe_names: "Probe names",
            sample_names: "Sample names",
        }
    }

    parameter_meta {
        idats: "Array of raw IDAT files from the Illumina methlyation array"
        disk_size_gb: "Disk size in GB"
    }

    input {
        Pair[File, File] idats
        Int disk_size_gb = 10
    }

    String idat_base = basename(idats.left)
    String out_base = sub(idat_base, "(_Grn|_Red).idat", "")

    #@ except: LineWidth
    command <<<
        echo "Processing IDAT files"
        ln -s ~{idats.left} ~{idats.right} .

        R --no-save <<SCRIPT
            options(error = function() traceback(3))

            library(minfi)
            library(dplyr)
            library(data.table)
            library(IlluminaHumanMethylationEPICmanifest)
            library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

            # These arrays have two types of probes (helpfully named Infinium type 1 and Infinium type 2).
            # The probes also cover a varying number of CpG sites.
            # The normalization method works by taking N probes of each type that each have 1, 2, or 3 CpGs.
            # So it selects 6N probes in total. This selection is "random" if the seed is not fixed, of course.
            # Since we're processing samples individually in this step instead of as a cohort,
            # we want the seed to be consistent so the same set of probes is chosen for each sample.
            set.seed(1)

            dir=getwd()

            RGSet <- read.metharray("~{idat_base}", verbose = TRUE, force = TRUE)
            saveRDS(RGSet, "~{out_base}.RGSet.rds")
            phenoData <- pData(RGSet)
            write.csv(phenoData, "~{out_base}.phenoData.csv")

            # The manifest is needed by preprocessRAW
            manifest <- getManifest(RGSet)
            manifest

            # Load raw data into a MethylSet object be converting red/green
            # channels into a matrix of methlyated and unmethylated signals.
            MSet <- preprocessRaw(RGSet)
            saveRDS(MSet, "~{out_base}.MSet.rds")

            # Convert to a RatioSet
            RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
            saveRDS(RSet, "~{out_base}.RSet.rds")

            # Add genomic coordinates to each probe (plus additional annotation)
            GRset <- mapToGenome(RSet)
            saveRDS(GRSet, "~{out_base}.GRSet.rds")

            # Take the genomic mapped RatioSet and fill Beta values (non-normalized).
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Get the NON-normalized beta values:
            beta <- getBeta(GRset)
            write.csv(beta, "~{out_base}.beta.csv")
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            # Get M-value matrix and copy-number matrix
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Get M and CN vals if interested:
            M <- getM(GRset)
            write.csv(M, "~{out_base}.m_values.csv")
            CN <- getCN(GRset)
            write.csv(CN, "~{out_base}.cn_values.csv")
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # Get sample names, probe names, and phenotype data
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            sampleNames <- sampleNames(GRset)
            write.csv(sampleNames, "~{out_base}.sampleNames.csv")
            probeNames <- featureNames(GRset)
            write.csv(probeNames, "~{out_base}.probeNames.csv")
            pheno <- pData(GRset)
            write.csv(pheno, "~{out_base}.pheno.csv")

            gr <- granges(GRset)
            write.csv(gr, "~{out_base}.gr.csv")

            annotation <- getAnnotation(GRset)
            write.csv(annotation, "~{out_base}.annotation.csv")
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Perform SWAN normalization on beta values
            GRset.swan_norm <- preprocessSWAN(RGSet)
            write.csv(GRset.swan_norm, "~{out_base}.GRset.swan_norm.csv")
            beta_swan_norm <- getBeta(GRset.swan_norm)
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Write the normalized beta-values that have NOT yet had
            # low-variance probes filtered out
            write.csv(beta_swan_norm, "~{out_base}.beta_swan_norm_unfiltered.csv")
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            # Write the normalized beta-values that have NOT yet had
            # low-variance probes filtered out
            # Filter to only those that are mappable to the genome.
            RSet <- ratioConvert(GRset.swan_norm)
            GRset <- mapToGenome(RSet)
            beta_swan_norm <- getBeta(GRset)
            write.csv(beta_swan_norm[order(rownames(beta_swan_norm)),], "~{out_base}.beta_swan_norm_unfiltered.genomic.csv")
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        SCRIPT
    >>>

    #@ except: LineWidth
    output {
        File beta_swan_norm_unfiltered = out_base + ".beta_swan_norm_unfiltered.csv"
        File beta_swan_norm_unfiltered_genomic = out_base + ".beta_swan_norm_unfiltered.genomic.csv"
        File mset = out_base + ".MSet.rds"
        File rgset = out_base + ".RGSet.rds"
        File rset = out_base + ".RSet.rds"
        File annotation = out_base + ".annotation.csv"
        File beta = out_base + ".beta.csv"
        File cn_values = out_base + ".cn_values.csv"
        File m_values = out_base + ".m_values.csv"
        File pheno = out_base + ".pheno.csv"
        File pheno_data = out_base + ".phenoData.csv"
        File probe_names = out_base + ".probeNames.csv"
        File sample_names = out_base + ".sampleNames.csv"
    }

    runtime {
        container: "ghcr.io/stjudecloud/minfi:branch-methylation-1.48.0-1"
        memory: "8 GB"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}
