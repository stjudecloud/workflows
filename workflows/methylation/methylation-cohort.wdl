version 1.1

workflow methylation_cohort {
    meta {
        description: "Process methylation data for a cohort of samples"
        warning: "We recommend against running this workflow direcly, and would suggest instead running the `methylation` workflow defined in `./methylation-standard.wdl`."
        outputs: {
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
        unfiltered_normalized_beta: "Array of unfiltered normalized beta values for each sample"
        p_values: "Array of detection p-value files for each sample."
        skip_pvalue_check: "Skip filtering based on p-values, even if `p_values` is supplied."
        num_probes: "Number of probes to use when filtering to the top `num_probes` probes with the highest standard deviation."
    }

    input {
        Array[File] unfiltered_normalized_beta
        Array[File] p_values = []
        Boolean skip_pvalue_check = false
        Int num_probes = 10000
    }

    Int max_length = 500
    Int beta_length = length(unfiltered_normalized_beta)
    Int pval_length = length(p_values)

    if (beta_length > max_length){
        scatter (merge_num in range((beta_length / max_length) + 1)){
            # Get the sublist of beta files
            scatter (beta_num in range(max_length)){
                Int num = (
                    if merge_num > 0
                    then beta_num + (merge_num * max_length)
                    else beta_num
                )
                if (num < beta_length){
                    File bam_list = unfiltered_normalized_beta[num]
                }
            }
        }
        scatter (iter_index in range(length(bam_list))){
            call combine_data as inner_merge { input:
                files_to_combine = select_all(bam_list[iter_index]),
                combined_file_name = "~{iter_index}.combined.csv",
            }
        }

        call combine_data as final_merge { input:
            files_to_combine = inner_merge.combined_file,
            combined_file_name = "combined_beta.csv",
        }

        if (pval_length > 0 && !skip_pvalue_check){
            # If p-values are provided, merge those as well
            scatter (merge_num in range((pval_length / max_length) + 1)){
                # Get the sublist of p-value files
                scatter (pval_num in range(max_length)){
                    Int num_p = (
                        if merge_num > 0
                        then pval_num + (merge_num * max_length)
                        else pval_num
                    )
                    if (num_p < pval_length){
                        File pval_list = p_values[num_p]
                    }
                }
            }
            scatter (iter_index in range(length(pval_list))){
                call combine_data as inner_merge_pvals { input:
                    files_to_combine = select_all(pval_list[iter_index]),
                    combined_file_name = "~{iter_index}.pvals.combined.csv",
                }
            }

            call combine_data as final_merge_pvals { input:
                files_to_combine = inner_merge_pvals.combined_file,
                combined_file_name = "combined_pvals.csv",
            }
        }
    }

    if (beta_length <= max_length){
        call combine_data as simple_merge { input:
            files_to_combine = unfiltered_normalized_beta,
            combined_file_name = "combined_beta.csv",
        }
        if (pval_length > 0 && !skip_pvalue_check){
            call combine_data as simple_merge_pval { input:
                files_to_combine = p_values,
                combined_file_name = "combined_pvals.csv",
            }
        }
    }

    File? pval_file = if !skip_pvalue_check then select_first(
        [
            final_merge_pvals.combined_file,
            simple_merge_pval.combined_file,
        ])
        else None

    call filter_probes { input:
        beta_values = select_first(
            [
                final_merge.combined_file,
                simple_merge.combined_file,
            ]
        ),
        p_values = pval_file,
        num_probes,
    }

    call generate_umap { input:
        filtered_beta_values = filter_probes.filtered_beta_values,
    }

    call plot_umap { input:
        umap = generate_umap.umap,
    }

    output {
        File combined_beta = select_first(
            [
                final_merge.combined_file,
                simple_merge.combined_file,
            ]
        )
        File filtered_beta = filter_probes.filtered_beta_values
        File filtered_probeset = filter_probes.filtered_probes
        File umap_embedding = generate_umap.umap
        File umap_plot = plot_umap.umap_plot
        File? probe_pvalues = pval_file
    }
}

task combine_data {
    meta {
        description: "Combine data from multiple CSV files by column"
        outputs: {
            combined_file: "Combined CSV file"
        }
    }

    parameter_meta {
        files_to_combine: "Array of files with values for each sample"
        combined_file_name: "Name of the combined file"
        simple_merge: "Use simple merge rather than batched read. Use this if different arrays are to be combined."
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] files_to_combine
        String combined_file_name = "combined.csv"
        Boolean simple_merge = false
        Int modify_memory_gb = 0
    }

    Int memory_gb = ceil(size(files_to_combine, "GiB") *
        if simple_merge then 2 else 1)
        + modify_memory_gb
        + 2
    Int disk_size_gb = ceil(size(files_to_combine, "GiB") * 2) + 2

    command <<<
        python /scripts/methylation/combine.py \
            --output-name "~{combined_file_name}" \
            ~{if simple_merge then "--simple-merge" else ""} \
            ~{sep(" ", quote(files_to_combine))}
    >>>

    output {
        File combined_file = combined_file_name
    }

    runtime {
        container: "ghcr.io/stjudecloud/pandas:2.2.1-6"
        memory: "~{memory_gb} GB"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}

task filter_probes {
    meta {
        description: "Filter probes based on standard deviation"
        help: "Probes are filtered by calculating the standard deviation of the beta value for each probe across all samples. The top `num_probes` probes with the highest standard deviation are retained."
        outputs: {
            filtered_beta_values: "Filtered beta values for all samples",
            filtered_probes: "Probes that were retained after filtering.",
        }
    }

    parameter_meta {
        beta_values: "Beta values for all samples"
        p_values: "P-values for all samples"
        prefix: "Prefix for the output files. The extensions `.beta.csv` and `.probes.csv` will be appended."
        pval_threshold: "P-value cutoff to determine poor quality probes"
        pval_sample_fraction: "Fraction of samples that must exceed p-value threshold to exclude probe"
        num_probes: "Number of probes to retain after filtering"
    }

    input {
        File beta_values
        File? p_values
        String prefix = "filtered"
        Float pval_threshold = 0.01
        Float pval_sample_fraction = 0.5
        Int num_probes = 10000
    }

    Int disk_size_gb = ceil(size(beta_values, "GiB") * 2) + 2

    command <<<
        python /scripts/methylation/filter.py \
            --output-name "~{prefix}.beta.csv" \
            --filtered-probes "~{prefix}.probes.csv" \
            --num-probes ~{num_probes} \
            --pval-threshold ~{pval_threshold} \
            --pval-sample-fraction ~{pval_sample_fraction} \
            ~{"--pval '" + p_values + "'"} \
            "~{beta_values}"
    >>>

    output {
        File filtered_beta_values = "~{prefix}.beta.csv"
        File filtered_probes = "~{prefix}.probes.csv"
    }

    runtime {
        container: "ghcr.io/stjudecloud/pandas:2.2.1-6"
        memory: "8 GB"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}

task generate_umap {
    meta {
        description: "Generate UMAP embedding"
        outputs: {
            umap: "UMAP embedding for all samples"
        }
    }

    parameter_meta {
        filtered_beta_values: "Filtered beta values for all samples"
        prefix: "Prefix for output file. The extension `.csv` will be appended."
    }

    input {
        File filtered_beta_values
        String prefix = "umap"
    }

    Int disk_size_gb = ceil(size(filtered_beta_values, "GiB") * 2) + 2

    command <<<
        python /scripts/methylation/generate_umap.py \
            --beta "~{filtered_beta_values}" \
            --output-name "~{prefix}.csv"
    >>>

    output {
        File umap = "~{prefix}.csv"
    }

    runtime {
        container: "ghcr.io/stjudecloud/umap:0.5.7-9"
        memory: "8 GB"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}

task plot_umap {
    meta {
        description: "Plot UMAP embedding"
        outputs: {
            umap_plot: "UMAP plot for all samples"
        }
    }

    parameter_meta {
        umap: "UMAP embedding for all samples"
        plot_file: "Name of the output plot file"
    }

    input {
        File umap
        String plot_file = "umap.png"
    }

    command <<<
        python /scripts/methylation/plot_umap.py \
            --umap "~{umap}" \
            --output-name "~{plot_file}"
    >>>

    output {
        File umap_plot = "~{plot_file}"
    }

    runtime {
        cpu: 1
        memory: "4 GB"
        disks: "4 GB"
        container: "ghcr.io/stjudecloud/python-plotting:2.0.5"
        maxRetries: 1
    }
}
