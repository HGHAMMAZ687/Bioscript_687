#!/bin/bash

# List of genes to filter
genes="nsp1 nsp2 nsp3 nsp4 nsp5 nsp6 nsp7 nsp8 nsp9 nsp10 nsp11 S N ns10"

# Input and output files
input_file="aa_frame_1.csv"
output_file="Non-synonymous_frame_1.csv"

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
