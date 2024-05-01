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
        restrictive: "If true, run a less permissive validation of field values. Otherwise, check against SAM spec-defined values."
    }

    input {
        ReadGroup read_group
        Array[String] required_fields = ["ID"]
        Boolean restrictive = true
    }

    # The SAM spec allows any printable ASCII character in header fields.
    String sam_spec_pattern = "[\\ -~]+"
    # We have the opinion that is too permissive for ID and SM.
    String id_pattern = if restrictive then "id" else sam_spec_pattern
    String sample_pattern = if restrictive then "sample.?" else sam_spec_pattern
    Array[String] platforms = ["CAPILLARY", "DNBSEQ", "ELEMENT", "HELICOS", "ILLUMINA", "IONTORRENT", "LS454", "ONT", "PACBIO", "SINGULAR", "SOLID", "ULTIMA"]

    command <<<
        error=0
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "ID") -eq 1 ]
        then
            if [ -z "~{read_group.ID}" ]
            then
                >&2 echo "ID is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.ID)}" == "true" ]
        then
            if [[ "~{read_group.ID}" =~ ^~{id_pattern}$ ]]
            then
                >&2 echo "ID must match pattern ~{id_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "SM") -eq 1 ]
        then
            if [ -z "~{read_group.SM}" ]
            then
                >&2 echo "SM is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.SM)}" == "true" ]
        then
            if [[ "~{read_group.SM}" =~ ^~{sample_pattern}$ ]]
            then
                >&2 echo "SM must match pattern ~{sample_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "BC") -eq 1 ]
        then
            if [ -z "~{read_group.BC}" ]
            then
                >&2 echo "BC is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.BC)}" == "true" ]
        then
            if [[ ! "~{read_group.BC}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "BC must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "CN") -eq 1 ]
        then
            if [ -z "~{read_group.CN}" ]
            then
                >&2 echo "CN is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.CN)}" == "true" ]
        then
            if [[ ! "~{read_group.CN}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "CN must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "DS") -eq 1 ]
        then
            if [ -z "~{read_group.DS}" ]
            then
                >&2 echo "DS is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.DS)}" == "true" ]
        then
            if [[ ! "~{read_group.DS}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "DS must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "DT") -eq 1 ]
        then
            if [ -z "~{read_group.DT}" ]
            then
                >&2 echo "DT is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.DT)}" == "true" ]
        then
            if [[ ! "~{read_group.DT}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "DT must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "FO") -eq 1 ]
        then
            if [ -z "~{read_group.FO}" ]
            then
                >&2 echo "FO is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.FO)}" == "true" ]
        then
            if [[ ! "~{read_group.FO}" =~ ^\*|[ACMGRSVTWYHKDBN]+$ ]]
            then
                >&2 echo "FO must match pattern \*|[ACMGRSVTWYHKDBN]+"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "KS") -eq 1 ]
        then
            if [ -z "~{if defined(read_group.KS) then read_group.KS else ""}" ]
            then
                >&2 echo "KS is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.KS)}" == "true" ]
        then
            if [[ ! "~{read_group.KS}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "KS must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "LB") -eq 1 ]
        then
            if [ -z "~{read_group.LB}" ]
            then
                >&2 echo "LB is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.LB)}" == "true" ]
        then
            if [[ ! "~{read_group.LB}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "LB must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PG") -eq 1 ]
        then
            if [ -z "~{read_group.PG}" ]
            then
                >&2 echo "PG is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.PG)}" == "true" ]
        then
            if [[ ! "~{read_group.PG}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "PG must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PI") -eq 1 ]
        then
            if [ -z "~{read_group.PI}" ]
            then
                >&2 echo "PI is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.PI)}" == "true" ]
        then
            if [[ ! "~{read_group.PI}" =~ ^[0-9]+$ ]]
            then
                >&2 echo "PI must match pattern [0-9]+"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PL") -eq 1 ]
        then
            if [ -z "~{read_group.PL}" ]
            then
                >&2 echo "PL is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.PL)}" == "true" ]
        then
            if [[ ! "~{read_group.PL}" =~ ^~{sep("|", platforms)}$ ]]
            then
                >&2 echo "PL must match pattern ~{sep("|", platforms)}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PM") -eq 1 ]
        then
            if [ -z "~{read_group.PM}" ]
            then
                >&2 echo "PM is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.PM)}" == "true" ]
        then
            if [[ ! "~{read_group.PM}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "PM must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi
        if [ $(echo "~{sep(" ", required_fields)}" | grep -Ewc "PU") -eq 1 ]
        then
            if [ -z "~{read_group.PU}" ]
            then
                >&2 echo "PU is required"
                error=1
            fi
        fi
        if [ "~{defined(read_group.PU)}" == "true" ]
        then
            if [[ ! "~{read_group.PU}" =~ ^~{sam_spec_pattern}$ ]]
            then
                >&2 echo "PU must match pattern ~{sam_spec_pattern}"
                error=1
            fi
        fi

        if [ $error -eq 1 ]
        then
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
        maxRetries: 0
    }
}