version 1.1

workflow methylation_cohort {
    meta {
        description: "Process methylation data for a cohort of samples"
        outputs: {
            combined_beta: "Matrix (in CSV format) containing beta values for every (common) probe on the array as rows and all of the input samples as columns.",
            filtered_beta: "Matrix (in CSV format) containing only beta values for the retained probes (top N highest standard deviation) for all provided samples.",
            filtered_probeset: "List of probe names that were retained after filtering to the top N highest standard deviation",
            umap_embedding: "UMAP embedding for all samples",
            umap_plot: "UMAP plot for all samples",
        }
    }

    parameter_meta {
        unfiltered_normalized_beta: "Array of unfiltered normalized beta values for each sample"
        num_probes: "Number of probes to use when filtering to the top `num_probes` probes with the highest standard deviation."
    }

    input {
        Array[File] unfiltered_normalized_beta
        Int num_probes = 10000
    }

    Int max_length = 500
    Int beta_length = length(unfiltered_normalized_beta)

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
                unfiltered_normalized_beta = select_all(bam_list[iter_index]),
                combined_file_name = "~{iter_index}.combined.csv",
                modify_memory_gb = 25,
            }
        }

        call combine_data as final_merge { input:
            unfiltered_normalized_beta = inner_merge.combined_beta,
        }
    }

    if (beta_length <= max_length){
        call combine_data as simple_merge { input:
            unfiltered_normalized_beta,
        }
    }

    call filter_probes { input:
        beta_values = select_first(
            [
                final_merge.combined_beta,
                simple_merge.combined_beta,
            ]
        ),
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
                final_merge.combined_beta,
                simple_merge.combined_beta,
            ]
        )
        File filtered_beta = filter_probes.filtered_beta_values
        File filtered_probeset = filter_probes.filtered_probes
        File umap_embedding = generate_umap.umap
        File umap_plot = plot_umap.umap_plot
    }
}

task combine_data {
    meta {
        description: "Combine data from multiple samples"
        outputs: {
            combined_beta: "Combined beta values for all samples"
        }
    }

    parameter_meta {
        unfiltered_normalized_beta: "Array of unfiltered normalized beta values for each sample"
        combined_file_name: "Name of the combined file"
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs. Specified in GB."
    }

    input {
        Array[File] unfiltered_normalized_beta
        String combined_file_name = "combined_beta.csv"
        Int modify_memory_gb = 0
    }

    Int memory_gb = ceil(size(unfiltered_normalized_beta, "GiB")) + modify_memory_gb + 2
    Int disk_size_gb = ceil(size(unfiltered_normalized_beta, "GiB") * 2) + 2

    command <<<
        python $(which combine.py) \
            --output-name ~{combined_file_name} \
            ~{sep(" ", unfiltered_normalized_beta)}
    >>>

    output {
        File combined_beta = combined_file_name
    }

    runtime {
        container: "ghcr.io/stjudecloud/pandas:branch-methylation-2.2.1-0"
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
        prefix: "Prefix for the output files. The extensions `.beta.csv` and `.probes.csv` will be appended."
        num_probes: "Number of probes to retain after filtering"
    }

    input {
        File beta_values
        String prefix = "filtered"
        Int num_probes = 10000
    }

    Int disk_size_gb = ceil(size(beta_values, "GiB") * 2) + 2

    command <<<
        python $(which filter.py) \
            --output-name ~{prefix}.beta.csv \
            --filtered-probes ~{prefix}.probes.csv \
            --num-probes ~{num_probes} \
            ~{beta_values}
    >>>

    output {
        File filtered_beta_values = "~{prefix}.beta.csv"
        File filtered_probes = "~{prefix}.probes.csv"
    }

    runtime {
        container: "ghcr.io/stjudecloud/pandas:branch-methylation-2.2.1-0"
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
        python $(which generate_umap.py) \
            --beta ~{filtered_beta_values} \
            --output-name ~{prefix}.csv
    >>>

    output {
        File umap = "~{prefix}.csv"
    }

    runtime {
        container: "ghcr.io/stjudecloud/umap:branch-methylation-0.5.7-0"
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
        python $(which plot_umap.py) --umap ~{umap} --output-name ~{plot_file}
    >>>

    output {
        File umap_plot = "~{plot_file}"
    }

    runtime {
        container: "ghcr.io/stjudecloud/python-plotting:branch-methylation-1.0.0"
        memory: "4 GB"
        cpu: 1
        disks: "4 GB"
        maxRetries: 1
    }
}
