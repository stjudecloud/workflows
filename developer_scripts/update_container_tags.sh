#!/usr/bin/env bash

set -euo pipefail

[[ "$#" -ne 1 ]] && >&2 echo "Usage: $(basename "$0") <target branch>" && exit 1
target=$1; shift;

if [[ "$target" == "main" ]]; then
  search_pattern='container: "ghcr\.io/stjudecloud/([^:]*):branch-[^-]*-(.*)"'
  replacement='container: "ghcr.io/stjudecloud/\1:\2"'
else
  search_pattern='container: "ghcr\.io/stjudecloud/([^:]*):(branch-[^-]*-)?(.*)"'
  replacement='container: "ghcr.io/stjudecloud/\1:branch-'$target'-\3"'
fi

# freeBSD sed (including the Mac version) requires a space and then an argument to -i
if [ "$(uname)" == "Darwin" ]
then
  inplace_arg='-i ""'
else
  inplace_arg='-i'
fi

find . -not -path '.git' -a -type f -name '*.wdl' | xargs sed -Er $inplace_arg "s,$search_pattern,$replacement,g"
