#!/bin/bash
# filepath: extract_params.sh

# ./extract_json.sh /path/to/json/files /path/to/output

# Check if arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    echo "Example: $0 ./input_json_dir ./output_dir"
    exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Set output file path
OUTPUT_FILE="$OUTPUT_DIR/results.csv"

# Create CSV file with headers
echo "BASE_NAME,PARAM_GR,RUN,NPAR,AIC,lnLik,PARAM_CONS,convergence" > "$OUTPUT_FILE"

# Find all JSON files in input directory and process them
find "$INPUT_DIR" -name "*.json" -exec jq -r '[
    .BASE_NAME[0],
    .PARAM_GR[0],
    .RUN[0],
    .NPAR[0],
    .AIC[0],
    .lnLik[0],
    .PARAM_CONS[0],
    .convergence[0]
] | @csv' {} \; >> "$OUTPUT_FILE"

echo "Extraction complete. Results saved to $OUTPUT_FILE"
