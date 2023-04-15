## # WDL tool template

version 1.0

task static_disk_and_ram_task {
    meta {
        description: "This template is appropriate for tasks with static disk space and RAM requirements. Appropriately update the default disk and RAM allocations." 
    }

    parameter_meta {
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        Int memory_gb = <>
        Int disk_size_gb = <>
        Int max_retries = 1
    }

    command <<<

    >>>

    output {

    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: ""
        maxRetries: max_retries
    }
}

task dynamic_disk_and_ram_task {
    meta {
        description: "This template is appropriate for tasks with dynamic disk and RAM requirements. Ensure the dynamic allocation of disk space and memory is sane." 
    }

    parameter_meta {
        modify_memory_gb: "Add to or subtract from dynamic memory allocation. Default memory is determined by the size of the inputs."
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        Int modify_memory_gb = 0
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Int input_size_gb = ceil(size(<input files>, "GiB"))

    Int disk_size_gb = ceil(input_size_gb * 1.5) + modify_disk_size_gb
    Int memory_gb = ceil(input_size_gb * 1.2) + modify_memory_gb

    command <<<

    >>>

    output {

    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: ""
        maxRetries: max_retries
    }
}

task detect_nproc_task {
    meta {
        description: "This template is appropriate for all tasks which can be run on multiple cores. Update the default disk and RAM allocations, or copy and paste from `dynamic_disk_and_ram_task` as appropriate." 
    }

    parameter_meta {
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task"
        ncpu: "Number of cores to allocate for task"
        detect_nproc: "Use all available cores? Recommended for cloud environments. Not recommended for cluster environments."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        Int memory_gb = <>
        Int disk_size_gb = <>
        Int ncpu = 1
        Boolean detect_nproc = false
        Int max_retries = 1
    }

    command <<<
        set -euo pipefail

        n_cores=~{ncpu}
        if [ "~{detect_nproc}" = "true" ]; then
            n_cores=$(nproc)
        fi
    >>>

    output {

    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        cpu: ncpu
        docker: ""
        maxRetries: max_retries
    }
}

task outfile_name_task {
    meta {
        description: "This template is appropriate for tasks where naming of the output file doesn't effect downstream analysis. Update the `parameter_meta` for `outfile_name` with the type of file in question, but do not change the variable name of `outfile_name`. Change `<output name>` to something short but descriptive. Appropriately set the default disk and RAM allocations." 
    }

    parameter_meta {
        outfile_name: "Name for the <type of file> file"
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        String outfile_name = basename(<input file>, ".<extension>") + ".<new extension>"
        Int memory_gb = <>
        Int disk_size_gb = <>
        Int max_retries = 1
    }

    command <<<

    >>>

    output {
        File <output name> = outfile_name
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: ""
        maxRetries: max_retries
    }
}

task prefix_task {
    meta {
        description: "This template is appropriate for tasks where the extension of the output file name effects downstream analysis. Update the `parameter_meta` for `prefix` with the type of file in question and the extension that will be added, but do not change the variable name of `prefix`. Change `<output name>` to something short but descriptive. Appropriately update the default disk and RAM allocations." 
    }

    parameter_meta {
        prefix: "Prefix for the <type of file> file. The extension `<extension>` will be added."
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        String prefix = basename(<input file>, ".<extension>") + ".<new extension>"
        Int memory_gb = 10
        Int disk_size_gb = 10
        Int max_retries = 1
    }

    command <<<

    >>>

    output {
        File <output name> = prefix + ".<extension>"
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: ""
        maxRetries: max_retries
    }
}

task string_choices_task {
    meta {
        description: "This template is appropriate for tasks that have a string parameter for which only >2 choices are valid. If there are 2 possible choices, consider constructing a Boolean instead. Task should fail quickly if an invalid choice is input." 
    }

    parameter_meta {
        <choice_input>: {
            description: "Description of the parameter"
            choices: [
                'foo',
                'bar',
                'baz'
            ]
        }
        memory_gb: "RAM to allocate for task"
        disk_size_gb: "Disk space to allocate for task"
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        Int memory_gb = <>
        Int disk_size_gb = <>
        Int max_retries = 1
    }

    command <<<

    >>>

    output {

    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size_gb + " GB"
        docker: ""
        maxRetries: max_retries
    }
}