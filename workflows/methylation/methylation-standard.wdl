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
            probes_with_snps: "List of probes that are affected by SNPs",
            sex_probe_list: "List of probes that map to sex chromosomes",
            high_pval_probes: "List of probes that were filtered out due to high p-values",
            non_genomic_probes: "List of probes that do not map to the genome",
        }
        allowNestedInputs: true
    }

    parameter_meta {
        green_idats: "Array of raw green IDAT files from the Illumina methylation array. See NOTE in `process_raw_idats` task for naming convention."
        red_idats: "Array of raw red IDAT files from the Illumina methylation array.  See NOTE in `process_raw_idats` task for naming convention."
        additional_probes_to_exclude: "Optional file containing a list of additional probes to exclude from the analysis. This should be a text file with one probe name per line."
    }

    input {
        Array[File] green_idats
        Array[File] red_idats
        File? additional_probes_to_exclude
    }

    scatter (pair in zip(green_idats, red_idats)) {
        call preprocess.process_raw_idats { input:
            idats = pair
        }
    }

    call preprocess.list_sex_probes {}

    call cohort.methylation_cohort { input:
        unfiltered_normalized_beta =
            process_raw_idats.beta_swan_norm_unfiltered_genomic,
        p_values = process_raw_idats.probe_pvalues,
        sex_probe_list = list_sex_probes.probe_list,
        additional_probes_to_exclude,
    }

    Array[File] probe_files = process_raw_idats.probes_with_snps
    Int probelist_length = length(probe_files)
    Int max_length = 100

    if (probelist_length > max_length){
        scatter (merge_num in range((probelist_length / max_length) + 1)){
            # Get the sublist of probe files
            scatter (probe_num in range(max_length)){
                Int num = (
                    if merge_num > 0
                    then probe_num + (merge_num * max_length)
                    else probe_num
                )
                if (num < probelist_length){
                    File probe_file_batches = probe_files[num]
                }
            }
        }
        scatter (iter_index in range(length(probe_file_batches))){
            call concat_and_uniq { input:
                files_to_combine = select_all(probe_file_batches[iter_index]),
                output_file_name = "probes_with_snps_part_~{iter_index}.tab",
            }
        }

        call concat_and_uniq as final_cat { input:
            files_to_combine = flatten([
                concat_and_uniq.combined_file
            ]),
            output_file_name = "probes_with_snps.tab",
        }
    }

    if (probelist_length <= max_length){
        call concat_and_uniq as simple_merge { input:
            files_to_combine = probe_files,
            output_file_name = "probes_with_snps.tab",
        }
    }

    Array[File] non_genomic_probe_list = process_raw_idats.non_genomic_probes
    Int non_genomic_probelist_length = length(non_genomic_probe_list)

    if (non_genomic_probelist_length > max_length){
        scatter (merge_num in range((non_genomic_probelist_length / max_length) + 1)){
            # Get the sublist of probe files
            scatter (probe_num in range(max_length)){
                Int num_ng = (
                    if merge_num > 0
                    then probe_num + (merge_num * max_length)
                    else probe_num
                )
                if (num_ng < non_genomic_probelist_length){
                    File non_genomic_probe_batches = non_genomic_probe_list[num_ng]
                }
            }
        }
        scatter (iter_index in range(length(non_genomic_probe_batches))){
            call concat_and_uniq as non_genomic_concat { input:
                files_to_combine = select_all(non_genomic_probe_batches[iter_index]),
                output_file_name = "non_genomic_probes_part_~{iter_index}.tab",
            }
        }

        call concat_and_uniq as final_cat_non_genomic { input:
            files_to_combine = flatten([
                non_genomic_concat.combined_file
            ]),
            output_file_name = "non_genomic_probes.tab",
        }
    }

    if (non_genomic_probelist_length <= max_length){
        call concat_and_uniq as simple_merge_non_genomic { input:
            files_to_combine = non_genomic_probe_list,
            output_file_name = "non_genomic_probes.tab",
        }
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
        File probes_with_snps = select_first([
            final_cat.combined_file,
            simple_merge.combined_file,
        ])
        File sex_probe_list = list_sex_probes.probe_list
        File? high_pval_probes = methylation_cohort.high_pval_probes
        File non_genomic_probes = select_first([
            final_cat_non_genomic.combined_file,
            simple_merge_non_genomic.combined_file,
        ])
    }
}

task concat_and_uniq {
    meta {
        description: "Concatenate multiple files and retain unique lines"
        outputs: {
            combined_file: "File containing unique lines from all input files"
        }
    }

    parameter_meta {
        files_to_combine: "Array of input files to concatenate"
        output_file_name: "Name of the output file"
    }

    input {
        Array[File] files_to_combine
        String output_file_name = "unique_lines.txt"
    }

    command <<<
        set -euo pipefail

        sort -u ~{sep(" ", quote(files_to_combine))} > "~{output_file_name}"
    >>>

    output {
        File combined_file = "~{output_file_name}"
    }

    runtime {
        container: "ghcr.io/stjudecloud/pandas:2.2.1-7"
        memory: "2 GB"
        cpu: 1
        disks: "4 GB"
        maxRetries: 1
    }
}
