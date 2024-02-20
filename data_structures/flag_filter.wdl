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

task validate_string_is_dec_oct_or_hex {
    meta {
        description: "Validates that a string is a decimal, octal, or hexadecimal number. **[WARNING]** Hexadecimal numbers must be prefixed with '0x' and only contain the characters [0-9A-F] to be valid (i.e. [a-f] is not allowed). Octal number must start with '0' and only contain the characters [0-7] to be valid. And decimal numbers must start with a digit between 1-9 and only contain the characters [0-9] to be valid."
        outputs: {
            is_valid: "True if the string is a decimal, octal, or hexadecimal number"
        }
    }

    parameter_meta {
        number: "The number to validate. See task description for accepted formats."
    }

    input {
        String number
    }

    command <<<
        if [[ "~{number}" =~ ^[1-9][0-9]+$ ]]; then
            # number is in decimal
            if [ "~{number}" -lt 4096 ]; then
                echo "Input number (~{number}) is valid" >> /dev/stderr
            else
                echo "Input number (~{number}) interpreted as decimal" >> /dev/stderr
                echo "But number must be less than 4096!" >> /dev/stderr
                exit 42
            fi
        elif [[ "~{number}" =~ ^0[0-7]{0,4}$ ]] \
            || [[ "~{number}" =~ ^0x[0-9A-F]{1,3}$ ]]
        then
            # number is in octal or hexadecimal
            # and number is less than 4096(decimal)
            echo "Input number (~{number}) is valid" >> /dev/stderr
        else
            # malformed for any reason
            echo "Input number (~{number}) is invalid" >> /dev/stderr
            echo "See task description for valid formats" >> /dev/stderr
            exit 42
        fi
    >>>

    output {
        String check = "passed"
    }

    runtime {
        memory: "4 GB"
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.3.0'
        maxRetries: 1
    }
}

workflow from_FlagFilterExplicit_to_FlagFilter {
    meta {
        description: "Converts a FlagFilterExplicit struct to a FlagFilter struct."
        output: {
            flag_filter: "FlagFilter struct"
        }
    }

    parameter_meta {
        flags: "FlagFilterExplicit struct to convert"
        validate_output: "If true, validate the output. This option is just for debugging purposes, and should be unnecassry in a production workflow."
    }

    input {
        FlagFilterExplicit flags
        Boolean validate_output = false
    }

    call from_BAMFlagsExplicit_to_String as include_if_any { input:
        flags = flags.include_if_any
    }
    call from_BAMFlagsExplicit_to_String as include_if_all { input:
        flags = flags.include_if_all
    }
    call from_BAMFlagsExplicit_to_String as exclude_if_any { input:
        flags = flags.exclude_if_any
    }
    call from_BAMFlagsExplicit_to_String as exclude_if_all { input:
        flags = flags.exclude_if_all
    }

    if (validate_output) {
        call validate_string_is_dec_oct_or_hex as validate_include_if_any { input:
            number = include_if_any.int_as_string
        }
        call validate_string_is_dec_oct_or_hex as validate_include_if_all { input:
            number = include_if_all.int_as_string
        }
        call validate_string_is_dec_oct_or_hex as validate_exclude_if_any { input:
            number = exclude_if_any.int_as_string
        }
        call validate_string_is_dec_oct_or_hex as validate_exclude_if_all { input:
            number = exclude_if_all.int_as_string
        }
    }

    FlagFilter result = FlagFilter {
        include_if_any: include_if_any.int_as_string,
        include_if_all: include_if_all.int_as_string,
        exclude_if_any: exclude_if_any.int_as_string,
        exclude_if_all: exclude_if_all.int_as_string,
    }

    output {
        FlagFilter flag_filter = result
    }
}
