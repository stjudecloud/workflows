version 1.1

import "./methylation-cohort.wdl" as cohort
import "./methylation-preprocess.wdl" as preprocess

workflow methylation {
    meta {
        name: "Methylation Standard"
        description: "Run the methylation pipeline on a cohort of samples"
        category: "Harmonization"
        outputs: {
            beta_swan_norm_unfiltered_genomic: "Normalized beta values for all probes that map to the genome",
            combined_beta: "Matrix (in CSV format) containing beta values for every (common) probe on the array as rows and all of the input samples as columns.",
            filtered_beta: "Matrix (in CSV format) containing only beta values for the retained probes (top N highest standard deviation) for all provided samples.",
            filtered_probeset: "List of probe names that were retained after filtering to the top N highest standard deviation",
            umap_embedding: "UMAP embedding for all samples",
            umap_plot: "UMAP plot for all samples",
            probe_pvalues: "Matrix (in CSV format) containing detection p-values for every (common) probe on the array as rows and all of the input samples as columns.",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        green_idats: "Array of raw green IDAT files from the Illumina methylation array. See NOTE in `process_raw_idats` task for naming convention."
        red_idats: "Array of raw red IDAT files from the Illumina methylation array.  See NOTE in `process_raw_idats` task for naming convention."
    }

    input {
        Array[File] green_idats
        Array[File] red_idats
    }

    scatter (pair in zip(green_idats, red_idats)) {
        call preprocess.process_raw_idats { input:
            idats = pair
        }
    }

    call cohort.methylation_cohort { input:
        unfiltered_normalized_beta =
            process_raw_idats.beta_swan_norm_unfiltered_genomic,
        p_values = process_raw_idats.probe_pvalues,
    }

    output {
        Array[File] beta_swan_norm_unfiltered_genomic =
            process_raw_idats.beta_swan_norm_unfiltered_genomic
        File combined_beta = methylation_cohort.combined_beta
        File filtered_beta = methylation_cohort.filtered_beta
        File filtered_probeset = methylation_cohort.filtered_probeset
        File umap_embedding = methylation_cohort.umap_embedding
        File umap_plot = methylation_cohort.umap_plot
        File? probe_pvalues = methylation_cohort.probe_pvalues
    }
}
