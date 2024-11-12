version 1.1

workflow methylation_cohort {
    meta {
        description: "Process methylation data for a cohort of samples"
        outputs: {
            combined_beta: "Combined beta values for all samples"
            filtered_beta: "Filtered beta values for all samples"
            umap_embedding: "UMAP embedding for all samples"
            umap_plot: "UMAP plot for all samples"
        }
    }

    parameter_meta {
        unfiltered_normalized_beta: "Array of unfiltered normalized beta values for each sample"
    }

    input {
        Array[File] unfiltered_normalized_beta
        Int max_length = 500
    }

    Int beta_length = length(unfiltered_normalized_beta)

    if (beta_length > max_length){
        scatter ( merge_num in range((beta_length / max_length) + 1)){
            # Get the sublist of beta files
            scatter ( beta_num in range(max_length)){
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
        scatter (iter in range(length(bam_list))){
            call combine_data as inner_merge { input:
                unfiltered_normalized_beta = select_all(bam_list[iter]),
                combined_file_name = iter + ".combined.csv"
            }
        }

        call combine_data as final_merge { input:
            unfiltered_normalized_beta = inner_merge.combined_beta,
        }
    }

    if (beta_length < max_length){
        call combine_data as simple_merge { input:
            unfiltered_normalized_beta,
        }
    }

    File merged_beta = select_first([final_merge.combined_beta, simple_merge.combined_beta])


    call filter_probes { input:
        combined_beta_values = merged_beta,
    }

    call generate_umap { input:
        filtered_beta_values = filter_probes.filtered_beta_values,
    }

    call plot_umap { input:
        umap = generate_umap.umap,
    }

    output {
        File combined_beta = merged_beta
        File filtered_beta = filter_probes.filtered_beta_values
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
    }

    input {
        Array[File] unfiltered_normalized_beta
        String combined_file_name = "combined_beta.csv"
    }

    Int disk_size_gb = ceil(size(unfiltered_normalized_beta, "GiB") * 2)

    command <<<
        echo "Combining data"
        input_files=""
        for file in ~{sep(" ", unfiltered_normalized_beta)}
        do
            #ln -s ~{sep(" ", unfiltered_normalized_beta)} .
            ln -s $file .
            input_files="${input_files} $(basename $file)"
        done

        cat <<SCRIPT > run.py
        import sys
        import argparse
        import pandas as pd

        def get_args():
            parser = argparse.ArgumentParser(
                description="Combine CSV files.")
            parser.add_argument(
                "csvs", type=str, nargs="+", help="List of CSV files.")

            args = parser.parse_args()

            return args
        
        def read(filename):
            df = pd.read_csv(filename)
            df.rename(columns={df.columns[0]: "probe"}, inplace=True)
            df.set_index("probe", inplace=True)
            return df

        if __name__ == "__main__":
            args = get_args()

            # Combine data
            df = pd.concat([read(f) for f in args.csvs], axis=1, join="inner")
            df.to_csv("~{combined_file_name}", index=True)
        SCRIPT

        python run.py $input_files

    >>>

    output {
        File combined_beta = combined_file_name
    }

    runtime {
        container: "quay.io/biocontainers/pandas:2.2.1"
        memory: "78 GB"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}

task filter_probes {
    meta {
        description: "Filter probes based on standard deviation"
        outputs: {
            filtered_beta_values: "Filtered beta values for all samples"
        }
    }

    parameter_meta {
        combined_beta_values: "Combined beta values for all samples"
    }

    input {
        File combined_beta_values
    }

    Int disk_size_gb = ceil(size(combined_beta_values, "GiB") * 2)

    command <<<
        ln -s ~{combined_beta_values} beta.csv

        python <<SCRIPT
        import pandas as pd

        # Read beta values
        beta = pd.read_csv("beta.csv", index_col=0)

        # Calculate standard deviation for each probe
        beta["sd"] = beta.std(axis='columns')

        # Filter to only the top 10,000 probes with the highest standard deviation
        beta.sort_values(by="sd", ascending=False).head(10000).drop(columns={"sd"}).to_csv("filtered_beta.csv")
        SCRIPT
    >>>

    output {
        File filtered_beta_values = "filtered_beta.csv"
    }

    runtime {
        container: "quay.io/biocontainers/pandas:2.2.1"
        memory: "178 GB"
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
    }

    input {
        File filtered_beta_values
    }

    Int disk_size_gb = ceil(size(filtered_beta_values, "GiB") * 2)

    command <<<
        ln -s ~{filtered_beta_values} beta.csv

        python <<SCRIPT
        import pandas as pd
        import umap

        # Read beta values
        beta = pd.read_csv("beta.csv", index_col=0)

        # Perform UMAP
        embedding = umap.UMAP().fit_transform(beta.T)

        # Save UMAP embedding
        umap = pd.DataFrame(data=embedding, index=beta.T.index, columns=["UMAP1", "UMAP2"])
        umap.index_name = "sample"
        umap.to_csv("umap.csv")
        SCRIPT
    >>>

    output {
        File umap = "umap.csv"
    }

    runtime {
        container: "ghcr.io/stjudecloud/umap:branch-methylation-0.5.7-1"
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
    }

    input {
        File umap
    }

    command <<<
        python <<SCRIPT
        import pandas as pd
        import matplotlib.pyplot as plt

        # Read UMAP embedding
        umap = pd.read_csv("umap.csv", index_col=0)

        # Plot UMAP
        plt.scatter(umap["UMAP1"], umap["UMAP2"])
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title("UMAP Embedding")
        plt.savefig("umap.png")    
        SCRIPT
    >>>

    output {
        File umap_plot = "umap.png"
    }

    runtime {
        container: "ghcr.io/stjudecloud/python-plotting:branch-methylation-1.0.0"
        memory: "4 GB"
        cpu: 1
        disks: "4 GB"
        maxRetries: 1
    }
}
