#!/usr/bin/env bash

set -euo pipefail

[[ "$#" -ne 1 ]] && >&2 echo "Usage: $(basename "$0") <target branch>" && exit 1
target=$1; shift;

if [[ "$target" == "main" ]]; then
  search_pattern='container: "ghcr\.io/stjudecloud/([^:]*):branch-[^-]*-(.*)"'
  replacement='container: "ghcr.io/stjudecloud/\1:\2"'
else
  search_pattern='container: "ghcr\.io/stjudecloud/([^:]*):(.*)"'
  replacement='container: "ghcr.io/stjudecloud/\1:branch-'$target'-\2"'
fi

find . -not -path '.git' -a -type f -name '*.wdl' | xargs sed -Eri '' "s,$search_pattern,$replacement,g"
