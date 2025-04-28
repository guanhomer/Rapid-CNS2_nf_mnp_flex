#!/bin/bash

# Input arguments
IN_FILE=$1
MNP_BED=$2
OUT_PATH=$3
FILENAME=$4  # id as filename

# Ensure output path exists
mkdir -p "$OUT_PATH"

# Filter for m (5mC) rows from bedMethyl file to prevent duplicated rows
awk '$4 == "m"' "$IN_FILE" > "${OUT_PATH}/${FILENAME}.tmp.bed"

# Intersect with the reference file for IlmnID using bedtools
bedtools intersect -a "${OUT_PATH}/${FILENAME}.tmp.bed" -b "$MNP_BED" -wa -wb > "${OUT_PATH}/${FILENAME}.tmp.bed"

# Add column names to the output file
column_names="chr start end coverage modC_percent IlmnID"
echo -e "$column_names" > "${OUT_PATH}/${FILENAME}.MNPFlex.subset.bed"

# Group by IlmnID (column $25) and summarize (sum) columns N_valid-cov ($10), N_mod ($12), and N_other-mod ($14)
# Methylation rate = ( N_mod + N_other-mod ) / N_valid-cov
awk -v FS="\t" -v OFS=" " '{
    coverage[$25] += $10; 
    modC[$25] += $12 + $14; 
    chr[$25] = $19; 
    start[$25] = $20; 
    end[$25] = $21
} 
END {
    for (id in coverage) {
        score = modC[id] / coverage[id] * 100
        # Print the result with space delimiters and 2 decimal places for the score
        printf "%s %s %s %d %.2f %s\n", chr[id], start[id], end[id], coverage[id], score, id
    }
}' "${OUT_PATH}/${FILENAME}.tmp.bed" >> "${OUT_PATH}/${FILENAME}.MNPFlex.subset.bed"

rm -r ${OUT_PATH}/${filename}.tmp.bed
