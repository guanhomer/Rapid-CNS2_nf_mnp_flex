#!/bin/bash

# Input arguments
IN_FILE=$1
MNP_BED=$2
OUT_PATH=$3
ID=$4  # ID as output filename

# Ensure output path exists
mkdir -p "$OUT_PATH"

# Filter for m (5mC) rows from bedMethyl file to prevent duplicated rows
# Intersect with the reference file for IlmnID using bedtools
awk '$4 == "m"' "$IN_FILE" | \
bedtools intersect -a stdin -b "$MNP_BED" -wa -wb > "${OUT_PATH}/${ID}.tmp.bed"

# Group by IlmnID (column $25) and summarize (sum) columns N_valid-cov ($10), N_mod ($12), and N_other-mod ($14)
# Methylation rate = ( N_mod + N_other-mod ) / N_valid-cov
# Create the output file with header and appended summary in one go
# Give pseudo count (-1234) for at least two MGMT probes
{
    echo "chr start end coverage methylation_percentage IlmnID"
    awk -v FS="\t" -v OFS=" " '
    {
        if ($25 == "MGMT") {
            mgmt_count++
            score = ($12 + $14) / $10 * 100
            printf "%s %s %s %d %.2f %s\n", $1, $2, $3, $10, score, $25
        } else {
            coverage[$25] += $10
            modC[$25] += $12 + $14
            chr[$25] = $19
            start[$25] = $20
            end[$25] = $21
        }
    }
    END {
        for (id in coverage) {
            score = modC[id] / coverage[id] * 100
            printf "%s %s %s %d %.2f %s\n", chr[id], start[id], end[id], coverage[id], score, id
        }
        if (mgmt_count < 2) {
            for (i = mgmt_count; i < 2; i++) {
                print "chr10 129467242 129467243 1 -0.001234 MGMT"
            }
        }
    }' "${OUT_PATH}/${ID}.tmp.bed"
} > "${OUT_PATH}/${ID}.MNPFlex.subset.bed"

rm ${OUT_PATH}/${ID}.tmp.bed
