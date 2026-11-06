#!/bin/bash

# ------- UNPARALLEL VERSION
# # Go to the directory containing the input files
# cd input

# # Loop through all .ini files and run treepl
# for file in input_*.ini
# do
#     echo "Processing $file..."
#     treepl "$file"
# done

# ------- PARALLEL VERSION
# Go to the directory containing the input files
mkdir -p timetrees

cd input

# Get number of CPU cores
num_cores=$(sysctl -n hw.ncpu)

# Use parallel to run TreePL jobs
# -j controls number of parallel jobs (using 75% of available cores)
ls input_*.ini | parallel -j$(($num_cores * 3/4)) 'echo "Processing {}..."; treepl {}'


