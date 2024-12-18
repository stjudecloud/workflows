version 1.1

workflow methylation_preprocess {
    meta {
        description: "Process raw IDAT files from the Illumina methylation array for a single sample"
        outputs: {
            beta_swan_norm_unfiltered: "Normalized beta values for all probes",
            beta_swan_norm_unfiltered_genomic: "Normalized beta values for all probes that map to the genome",
            annotation: "Annotation table for probes. Contains genomic coordinates and sequence and other information about the probes.",
            beta_unnorm: "Non-normalized beta values",
            cn_values: "Copy number values",
            m_values: "M values",
            probe_names: "Probe names found on the array",
            sample_names: "Sample name found on the array",
        }
    }

    parameter_meta {
        idats: "Pair of raw IDAT files from the Illumina methylation array"
        seed: "Random number generator seed for reproducibility. If processing samples independently, this should remain fixed for all runs that will be compared as a cohort."
    }

    input {
        Pair[File, File] idats
        Int seed = 1
    }

    call process_raw_idats { input:
        idats,
        seed,
    }

    #@ except: LineWidth
    output {
        File beta_swan_norm_unfiltered
            = process_raw_idats.beta_swan_norm_unfiltered
        File beta_swan_norm_unfiltered_genomic
            = process_raw_idats.beta_swan_norm_unfiltered_genomic
        File annotation = process_raw_idats.annotation
        File beta_unnorm = process_raw_idats.beta_unnorm
        File cn_values = process_raw_idats.cn_values
        File m_values = process_raw_idats.m_values
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
            annotation: "Annotation table for probes. Contains genomic coordinates and sequence and other information about the probes.",
            beta_unnorm: "Non-normalized beta values",
            cn_values: "Copy number values",
            m_values: "M values",
            probe_names: "Probe names found on the array",
            sample_names: "Sample name found on the array",
        }
    }

    parameter_meta {
        idats: "Array of raw IDAT files from the Illumina methlyation array"
        seed: "Random number generator seed for reproducibility. If processing samples independently, this should remain fixed for all runs that will be compared as a cohort."
        disk_size_gb: "Disk size in GB"
    }

    input {
        Pair[File, File] idats
        Int seed = 1
        Int disk_size_gb = 10
    }

    String idat_base = basename(idats.left)
    String out_base = sub(idat_base, "(_Grn|_Red).idat", "")

    command <<<
        set -euo pipefail

        ln -s ~{idats.left} ~{idats.right} .

        Rscript $(which methylation-preprocess.R) \
            --idat_base ~{idat_base} \
            --out_base ~{out_base} \
            --seed ~{seed}
    >>>

    output {
        File beta_swan_norm_unfiltered
            = out_base + ".beta_swan_norm_unfiltered.csv"
        File beta_swan_norm_unfiltered_genomic
            = out_base + ".beta_swan_norm_unfiltered.genomic.csv"
        File annotation = out_base + ".annotation.csv"
        File beta_unnorm = out_base + ".beta.csv"
        File cn_values = out_base + ".cn_values.csv"
        File m_values = out_base + ".m_values.csv"
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
