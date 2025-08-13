#!/usr/bin/env bash

set -euo pipefail

usage_and_quit () {
    echo "Runs sprocket or miniwdl (as determined by the RUNNER env var) with the supplied arguments"
    echo
    echo "If RUNNER env var is unset, sprocket will be run"
    echo
    echo "Usage: RUNNER=<sprocket or miniwdl> $(basename "$0") [path to WDL document] [input files and key value pairs]..." 
    exit 1
}

[[ "$#" -eq 0 ]] && usage_and_quit
[[ "$#" -le 2 ]] && usage_and_quit

runner=${RUNNER:-sprocket}
if [ "$runner" != "sprocket" ] && [ "$runner" != "miniwdl" ]; then
    echo "RUNNER environment variable must be 'sprocket' or 'miniwdl'"
    echo
    usage_and_quit
fi

wdl=$1; shift;
