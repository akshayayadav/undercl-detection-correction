#!/bin/bash


: ${num_cores:=1}

while getopts ":c:" opt; do
  case $opt in
    c) num_cores="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

snakemake -s /usr/local/bin/snakefile_subset-F1 -j $num_cores

