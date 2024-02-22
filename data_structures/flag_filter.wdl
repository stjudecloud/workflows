version 1.1

struct FlagFilter {
    String include_if_any
    String include_if_all
    String exclude_if_any
    String exclude_if_all
}

task validate_string_is_oct_dec_or_hex {
    meta {
        description: "Validates that a string is a decimal, octal, or hexadecimal number and less than 2^12. **[WARNING]** Hexadecimal numbers must be prefixed with '0x' and only contain the characters [0-9A-F] to be valid (i.e. [a-f] is not allowed). Octal number must start with '0' and only contain the characters [0-7] to be valid. And decimal numbers must start with a digit between 1-9 and only contain the characters [0-9] to be valid."
        outputs: {
            check: "Dummy output to enable caching."
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

workflow validate_FlagFilter {
    meta {
        description: "Validates a FlagFilter struct."
        output: {
            check: "Dummy output to enable caching."
        }
    }

    parameter_meta {
        flags: "FlagFilter struct to validate"
    }

    input {
        FlagFilter flags
    }

    call validate_string_is_oct_dec_or_hex as validate_include_if_any { input:
        number = flags.include_if_any
    }
    call validate_string_is_oct_dec_or_hex as validate_include_if_all { input:
        number = flags.include_if_all
    }
    call validate_string_is_oct_dec_or_hex as validate_exclude_if_any { input:
        number = flags.exclude_if_any
    }
    call validate_string_is_oct_dec_or_hex as validate_exclude_if_all { input:
        number = flags.exclude_if_all
    }

    output {
        String check = "passed"
    }

}
