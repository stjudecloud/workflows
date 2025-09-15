version 1.1

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
            probe_pvalues: "Matrix (in CSV format) containing detection p-values for every (common) probe on the array as rows.",
        }
    }

    parameter_meta {
        idats: "Array of raw IDAT files from the Illumina methlyation array. NOTE: these files must follow the normal naming convention for Illumina array data, that is the sample name followed by the green or red channel (e.g. `5723646052_R02C02_Grn.idat` and `5723646052_R02C02_Red.idat`)."
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

        ln -s "~{idats.left}" "~{idats.right}" .

        Rscript /scripts/methylation/methylation-preprocess.R \
            --idat_base "~{idat_base}" \
            --out_base "~{out_base}" \
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
        File probe_pvalues = out_base + ".detectionP.csv"
    }

    runtime {
        container: "ghcr.io/stjudecloud/minfi:branch-python-1.48.0-7"
        memory: "8 GB"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}
