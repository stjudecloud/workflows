#!/bin/bash

shopt -s extglob

for dir in docker/!(cellranger); do
  echo $dir
  tool_name=$(basename "$(echo "$dir" | awk '{print tolower($0)}')")
  echo $tool_name
  for tag in "$dir"/*; do
    tag=$(basename "$tag")
    echo $tag
  done
done

