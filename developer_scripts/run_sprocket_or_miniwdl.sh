#!/usr/bin/env bash

set -euo pipefail

usage_and_quit () {
    echo "Runs sprocket or miniwdl (as determined by the RUNNER env var) with the supplied arguments"
    echo
    echo "If RUNNER env var is unset, sprocket will be run"
    echo
    echo "Usage: RUNNER=<sprocket or miniwdl> $(basename "$0") [-i <input file>] [-t|-w <task or workflow (respectively) to run>] [path to WDL document] [optional key value pairs]..." 
    exit 1
}

[[ "$#" -eq 0 ]] && usage_and_quit

runner=${RUNNER:-sprocket}
if [ "$runner" != "sprocket" ] && [ "$runner" != "miniwdl" ]; then
    echo "ERROR: RUNNER environment variable must be 'sprocket' or 'miniwdl'"
    echo
    usage_and_quit
fi

input_file=""
task=""
wf=""
while getopts "hi:t:w:" option; do
   case $option in
      i)
         input_file=$OPTARG
         ;;
      t)
         task=$OPTARG
         ;;
      w)
         wf=$OPTARG
         ;;
      h)
         usage_and_quit
         ;;
      *)
         echo "ERROR: unrecognized flag"
         echo
         usage_and_quit
         ;;
   esac
done

shift $((OPTIND-1))

wdl=$1; shift;

if [[ $task ]] && [[ $wf ]]; then
   echo "ERROR: both -t and -w cannot be supplied"
   echo
   usage_and_quit
fi
if [ -z "$task" ] && [ -z "$wf" ]; then
   echo "ERROR: -t or -w must be supplied"
   echo
   usage_and_quit
fi

if [ "$runner" != "miniwdl" ]; then
   if [[ $task ]]; then
      entrypoint=$task
   else
      entrypoint=$wf
   fi
   sprocket run --output output --overwrite -e "$entrypoint" "$wdl" ${input_file:+"$input_file"} "$@"
else
   if [[ $input_file ]]; then
      input_dir=$(dirname "$input_file")
      sed -e 's;../../;'"$input_dir"'/../../;g' "$input_file" > tmp.json
      input_file=tmp.json
   fi
   if [ -d output ]; then
      rm -rf output
   fi
   miniwdl run ${task:+--task "$task"} --verbose --dir output/. ${input_file:+-i "$input_file"} "$wdl" "$@"
   if [ -e tmp.json ]; then
      rm tmp.json
   fi
fi
