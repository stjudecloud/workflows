version 1.1

# See the `read_groups` `parameter_meta` for definitions of each field
struct ReadGroup {
    String ID
    String? BC
    String? CN
    String? DS
    String? DT
    String? FO
    String? KS
    String? LB
    String? PG
    Int? PI
    String? PL
    String? PM
    String? PU
    String? SM
}

task ReadGroup_to_string {
    meta {
        description: "Stringifies a ReadGroup struct"
        outputs: {
            stringified_read_group: "Input ReadGroup as a string"
        }
    }

    parameter_meta {
        read_group: "ReadGroup struct to stringify"
    }

    input {
        ReadGroup read_group
    }

    command <<<
        {
            # TODO: I think this can be simplified by dropping the `if defined` checks?
            echo -n "~{'ID:~{read_group.ID}'}"  # required field. All others optional
            echo -n "~{if defined(read_group.BC) then ' BC:~{read_group.BC}' else ''}"
            echo -n "~{if defined(read_group.CN) then ' CN:~{read_group.CN}' else ''}"
            echo -n "~{if defined(read_group.DS) then ' DS:~{read_group.DS}' else ''}"
            echo -n "~{if defined(read_group.DT) then ' DT:~{read_group.DT}' else ''}"
            echo -n "~{if defined(read_group.FO) then ' FO:~{read_group.FO}' else ''}"
            echo -n "~{if defined(read_group.KS) then ' KS:~{read_group.KS}' else ''}"
            echo -n "~{if defined(read_group.LB) then ' LB:~{read_group.LB}' else ''}"
            echo -n "~{if defined(read_group.PG) then ' PG:~{read_group.PG}' else ''}"
            echo -n "~{if defined(read_group.PI) then ' PI:~{read_group.PI}' else ''}"
            echo -n "~{if defined(read_group.PL) then ' PL:~{read_group.PL}' else ''}"
            echo -n "~{if defined(read_group.PM) then ' PM:~{read_group.PM}' else ''}"
            echo -n "~{if defined(read_group.PU) then ' PU:~{read_group.PU}' else ''}"
            echo "~{if defined(read_group.SM) then ' SM:~{read_group.SM}' else ''}"
        } > out.txt
    >>>

    output {
        String stringified_read_group = read_string("out.txt")
    }

    runtime {
        memory: "4 GB"
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.4.0'
        maxRetries: 1
    }
}

task get_ReadGroups {
    meta {
        description: "Gets read group information from a BAM file and writes it out as JSON which is converted to a WDL struct."
        outputs: {
            read_groups: "An array of strings containing read group information. If `format_for_star = true`, all found read groups are contained in one string (`read_groups[0]`). If `format_for_star = false`, each found @RG line will be its own entry in output array `read_groups`."
        }
    }

    parameter_meta {
        bam: "Input BAM format file to get read groups from"
        modify_disk_size_gb: "Add to or subtract from dynamic disk space allocation. Default disk size is determined by the size of the inputs. Specified in GB."
    }

    input {
        File bam
        Int modify_disk_size_gb = 0
    }

    Float bam_size = size(bam, "GiB")
    Int disk_size_gb = ceil(bam_size) + 10 + modify_disk_size_gb

    command <<<
        set -euo pipefail
        BAM="~{bam}" OUTFILE="read_groups.json" python - <<END
import os  # lint-check: ignore
import pysam  # lint-check: ignore
import json  # lint-check: ignore
sam = pysam.AlignmentFile(os.environ["BAM"], "rb")

out_file = open(os.environ["OUTFILE"], "w")
json.dump(sam.header.to_dict()["RG"], out_file)
out_file.close()
END
    >>>

    output {
        Array[ReadGroup] read_groups = read_json("read_groups.json")
    }

    runtime {
        memory: "4 GB"
        disk: "~{disk_size_gb} GB"
        container: 'quay.io/biocontainers/pysam:0.22.0--py38h15b938a_1'
        maxRetries: 1
    }
}

task validate_ReadGroup {
    meta {
        description: "Validate a ReadGroup struct's fields are defined"
        outputs: {
            check: "Dummy output to indicate success and enable call-caching"
        }
    }

    parameter_meta {
        read_group: "ReadGroup struct to validate"
        required_fields: "Array of read group fields that must be defined"
    }

    input {
        ReadGroup read_group
        Array[String] required_fields = ["ID"]
    }

    command <<<
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "ID") -eq 1 ] && [ -z "~{read_group.ID}" ]; then
            >&2 echo "ID is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "SM") -eq 1 ] && [ -z "~{read_group.SM}" ]; then
            >&2 echo "SM is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "BC") -eq 1 ] && [ -z "~{read_group.BC}" ]; then
            >&2 echo "BC is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "CN") -eq 1 ] && [ -z "~{read_group.CN}" ]; then
            >&2 echo "CN is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "DS") -eq 1 ] && [ -z "~{read_group.DS}" ]; then
            >&2 echo "DS is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "DT") -eq 1 ] && [ -z "~{read_group.DT}" ]; then
            >&2 echo "DT is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "FO") -eq 1 ] && [ -z "~{read_group.FO}" ]; then
            >&2 echo "FO is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "KS") -eq 1 ] && [ -z "~{read_group.KS}" ]; then
            >&2 echo "KS is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "LB") -eq 1 ] && [ -z "~{read_group.LB}" ]; then
            >&2 echo "LB is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PG") -eq 1 ] && [ -z "~{read_group.PG}" ]; then
            >&2 echo "PG is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PI") -eq 1 ] && [ -z "~{read_group.PI}" ]; then
            >&2 echo "PI is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PL") -eq 1 ] && [ -z "~{read_group.PL}" ]; then
            >&2 echo "PL is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PM") -eq 1 ] && [ -z "~{read_group.PM}" ]; then
            >&2 echo "PM is required"
            exit 1
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PU") -eq 1 ] && [ -z "~{read_group.PU}" ]; then
            >&2 echo "PU is required"
            exit 1
        fi
    >>>

    output {
        String check = "passed"
    }

    runtime {
        memory: "4 GB"
        disk: "10 GB"
        container: 'ghcr.io/stjudecloud/util:1.4.0'
        maxRetries: 1
    }
}