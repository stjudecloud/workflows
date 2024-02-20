version 1.1

struct BAMFlagsExplicit {
    Boolean segmented  # 0x1
    Boolean segments_properly_aligned  # 0x2
    Boolean unmapped  # 0x4
    Boolean mate_unmapped  # 0x8
    Boolean reverse_complimented  # 0x10
    Boolean mate_reverse_complimented  # 0x20
    Boolean first  # 0x40
    Boolean last  # 0x80
    Boolean secondary  # 0x100
    Boolean qcfail  # 0x200
    Boolean duplicate  # 0x400
    Boolean supplementary  # 0x800
}

struct FlagFilterExplicit {
    BAMFlagsExplicit include_if_any
    BAMFlagsExplicit include_if_all
    BAMFlagsExplicit exclude_if_any
    BAMFlagsExplicit exclude_if_all
}

struct FlagFilter {
    String include_if_any
    String include_if_all
    String exclude_if_any
    String exclude_if_all
}

task from_BAMFlagsExplicit_to_String {
    meta {
        description: "Converts a BAMFlagsExplicit struct to an integer and stores it in a string"
        outputs: {
            int_as_string: "Input BAMFlagsExplicit as a string"
        }
    }

    parameter_meta {
        flags: "BAMFlagsExplicit struct to stringify"
    }

    input {
        BAMFlagsExplicit flags
    }

    command <<<
        BINARY_NUM=$(
            echo -n "~{if flags.supplementary then 1 else 0}"  # 0x800
            echo -n "~{if flags.duplicate then 1 else 0}"  # 0x400
            echo -n "~{if flags.qcfail then 1 else 0}"  # 0x200
            echo -n "~{if flags.secondary then 1 else 0}"  # 0x100
            echo -n "~{if flags.last then 1 else 0}"  # 0x80
            echo -n "~{if flags.first then 1 else 0}"  # 0x40
            echo -n "~{if flags.mate_reverse_complimented then 1 else 0}"  # 0x20
            echo -n "~{if flags.reverse_complimented then 1 else 0}"  # 0x10
            echo -n "~{if flags.mate_unmapped then 1 else 0}"  # 0x8
            echo -n "~{if flags.unmapped then 1 else 0}"  # 0x4
            echo -n "~{if flags.segments_properly_aligned then 1 else 0}"  # 0x2
            echo "~{if flags.segmented then 1 else 0}"  # 0x1
        )
        echo "$((2#$BINARY_NUM))" > out.txt
    >>>

    output {
        String int_as_string = read_string("out.txt")
    }

    runtime {
        memory: "4 GB"
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}

workflow from_FlagFilterExplicit_to_FlagFilter {
    input {
        FlagFilterExplicit flags
    }

    call from_BAMFlagsExplicit_to_String as include_if_any { input: flags = flags.include_if_any }
    call from_BAMFlagsExplicit_to_String as include_if_all { input: flags = flags.include_if_all }
    call from_BAMFlagsExplicit_to_String as exclude_if_any { input: flags = flags.exclude_if_any }
    call from_BAMFlagsExplicit_to_String as exclude_if_all { input: flags = flags.exclude_if_all }

    FlagFilter result = FlagFilter {
        include_if_any: include_if_any.int_as_string,
        include_if_all: include_if_all.int_as_string,
        exclude_if_any: exclude_if_any.int_as_string,
        exclude_if_all: exclude_if_all.int_as_string
    }

    output {
        FlagFilter flag_filter = result
    }
}
