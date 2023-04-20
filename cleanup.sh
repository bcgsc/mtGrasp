#!/bin/bash
set -eu -o pipefail


if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <DIRECTORY>"
  exit 1
fi


dir="$1"
cd "$dir"

# Delete the directories, printing out the full paths of the deleted directories
for d in abyss assembly_output blast bloom_filter_calc end_recovery mito_filtering_output pilon prepolishing subsets; do
  if [ -d "$d" ]; then
    path=$(realpath "$d")
    rm -r "$d" && echo "Removed directory: $path"
  fi
done