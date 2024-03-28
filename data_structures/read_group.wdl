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
        container: 'ghcr.io/stjudecloud/util:1.3.0'
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
