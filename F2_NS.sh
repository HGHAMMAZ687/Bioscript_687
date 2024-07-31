#!/bin/bash

# List of genes to filter
genes="M ns7b ns8"

# Input and output files
input_file="aa_frame_2.csv"
output_file="Non-synonymous_frame_2.csv"

# Run awk to filter the CSV
awk -F, -v genes="$genes" '
BEGIN {
    split(genes, gene_array, " ")
    for (i in gene_array) {
        gene_map[gene_array[i]] = 1
    }
}
NR == 1 || ($9 == "Non-synonymous" && $10 in gene_map && $8 != "?")
' "$input_file" > "$output_file"

echo "Filtered mutations have been written to $output_file"