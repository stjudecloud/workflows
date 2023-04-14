## WDL tool template

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
        Int memory_gb = 10
        Int disk_size_gb = 10
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

task dynamic_disk_task {
    meta {
        description: "This template is appropriate for tasks with large inputs and a static RAM requirement. Ensure the dynamic allocation of disk space is sane and appropriately update the default RAM allocation." 
    }

    parameter_meta {
        memory_gb: "RAM to allocate for task"
        modify_disk_size_gb: "Add or subtract disk space from dynamic allocation. Default disk size is determined by the size of the inputs."
        max_retries: "Number of times to retry in case of failure"
    }

    input {
        Int memory_gb = 10
        Int modify_disk_size_gb = 0
        Int max_retries = 1
    }

    Int disk_size_gb = ceil(size(<input files>, "GiB") * 1.5) + modify_disk_size_gb

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
        modify_memory_gb: "Add or subtract memory from dynamic allocation. Default memory is determined by the size of the inputs."
        modify_disk_size_gb: "Add or subtract disk space from dynamic allocation. Default disk size is determined by the size of the inputs."
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