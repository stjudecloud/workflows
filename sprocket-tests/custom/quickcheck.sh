#!/bin/bash

set -euo pipefail

out_json=$1

out_bam=$(jq -r .bam "$out_json")

samtools quickcheck "$out_bam"

