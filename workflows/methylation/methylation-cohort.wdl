version 1.1

workflow methylation_cohort {
    meta {
        description: "Process methylation data for a cohort of samples"
        outputs: {
            combined_beta: "Matrix (in CSV format) containing beta values for every (common) probe on the array as rows and all of the input samples as columns.",
            filtered_beta: "Matrix (in CSV format) containing only beta values for the retained probes (top N highest standard deviation) for all provided samples.",
            filtered_probes: "List of probe names that were retained after filtering to the top N highest standard deviation",
            umap_embedding: "UMAP embedding for all samples",
            umap_plot: "UMAP plot for all samples",
        }
    }

    parameter_meta {
        unfiltered_normalized_beta: "Array of unfiltered normalized beta values for each sample"
        max_length: "Maximum number of beta files to merge before using iteration"
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
                combined_file_name = iter + ".combined.csv",
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
        File filtered_probes = filter_probes.filtered_probes
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

    Int memory_gb = ceil(size(unfiltered_normalized_beta, "GiB") * 2) + addl_memory_gb
    Int disk_size_gb = ceil(size(unfiltered_normalized_beta, "GiB") * 2)

    command <<<
        echo "Combining data"
        input_files=""
        for file in ~{sep(" ", unfiltered_normalized_beta)}
        do
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

    #@ except: ContainerValue
    runtime {
        container: "quay.io/biocontainers/pandas:2.2.1"
        memory: "~{memory_gb} GB"
        cpu: 1
        disks: "~{disk_size_gb} GB"
        maxRetries: 1
    }
}

task filter_probes {
    meta {
        description: "Filter probes based on standard deviation"
        outputs: {
            filtered_beta_values: "Filtered beta values for all samples",
            filtered_probes: "Probes that were retained after filtering.",
        }
    }

    parameter_meta {
        beta_values: "Beta values for all samples"
        num_probes: "Number of probes to retain after filtering"
    }

    input {
        File beta_values
        Int num_probes = 10000
    }

    Int disk_size_gb = ceil(size(beta_values, "GiB") * 2)

    #@ except: LineWidth
    command <<<
        ln -s ~{beta_values} beta.csv

        cat <<SCRIPT > filter.py
        import csv
        import pandas as pd
        import numpy as np
        import sys
        import argparse

        def get_args():
            parser = argparse.ArgumentParser(
                description="Filter probes based on standard deviation.")
            parser.add_argument(
                "--num_probes", type=int, default=10000, help="Number of probes to retain after filtering.")
            parser.add_argument(
                "beta", type=str, help="Beta values CSV file.")

            args = parser.parse_args()

            return args

        if __name__ == "__main__":
            args = get_args()

            # Read beta values and compute standard deviation
            data = []
            with open("beta.csv", "r") as f:
                reader = csv.reader(f)
                header = next(reader)
                for i, line in enumerate(reader):
                    probe = line[0]
                    sd = np.std([float(x) for x in line[1:]])
                    data.append([probe, sd])

            sd_df = pd.DataFrame(data, columns=["probe", "sd"]).set_index("probe")

            # Filter probes based on standard deviation
            filtered_probes = sd_df.sort_values("sd", ascending=False).head(args.num_probes).index

            pd.Series(filtered_probes, index=filtered_probes).to_csv('filtered_probes.csv', index=False)

            # Filter beta values
            with open("beta.csv", "r") as f:
                reader = csv.reader(f)
                header = next(reader)
                header[0] = "probe"
                beta_data = [line for line in reader if line[0] in filtered_probes]

            bd = pd.DataFrame(beta_data, columns=header).set_index("probe")
            bd.to_csv("filtered_beta.csv")
        SCRIPT

        python filter.py --num_probes ~{num_probes} beta.csv
    >>>

    output {
        File filtered_beta_values = "filtered_beta.csv"
        File filtered_probes = "filtered_probes.csv"
    }

    #@ except: ContainerValue
    runtime {
        container: "quay.io/biocontainers/pandas:2.2.1"
        memory: "28 GB"
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

    #@ except: LineWidth
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

    #@ except: ContainerValue
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
        umap = pd.read_csv("~{umap}", index_col=0)

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

    #@ except: ContainerValue
    runtime {
        container: "ghcr.io/stjudecloud/python-plotting:branch-methylation-1.0.0"
        memory: "4 GB"
        cpu: 1
        disks: "4 GB"
        maxRetries: 1
    }
}
