#!/usr/bin/env bash

set -euo pipefail

[[ "$#" -ne 1 ]] && >&2 echo "Usage: $(basename "$0") <target branch>" && exit 1
target=$1; shift;

# Convert repository owner to lowercase for Docker image tags
repo_owner_lower=$(echo "${GITHUB_REPOSITORY_OWNER:-stjudecloud}" | tr '[:upper:]' '[:lower:]')

if [[ "$target" == "main" ]]; then
  search_pattern='container: "ghcr\.io/[^/]+/([^:]*):branch-[^-]*-(.*)"'
  replacement='container: "ghcr.io/'$repo_owner_lower'/\1:\2"'
else
  search_pattern='container: "ghcr\.io/[^/]+/([^:]*):(branch-[^-]*-)?(.*)"'
  replacement='container: "ghcr.io/'$repo_owner_lower'/\1:branch-'$target'-\3"'
fi

# freeBSD sed (including the Mac version) requires a space and then an argument to -i
if [ "$(uname)" == "Darwin" ]
then
  find . -not -path '.git' -a -type f -name '*.wdl' | xargs sed -Er -i '' "s,$search_pattern,$replacement,g"
else
  find . -not -path '.git' -a -type f -name '*.wdl' | xargs sed -Er -i'' "s,$search_pattern,$replacement,g"
fi
