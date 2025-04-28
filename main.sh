#!/bin/bash

set -euo pipefail

# Inputs
bam_dir="/path/to/folder/of/aligned/bams"             # Directory containing your BAM files
ref="/path/to/GRCh38.p14.genome.fa"
id="sample_name"
modkit_threads=8
WORKFLOW_ROOT="/path/to/workflow"

# Related data and scripts
MNPFLEX_BED="${WORKFLOW_ROOT}/data/MNP-flex.bed"
MNPFLEX_REPORT_SCRIPT="${WORKFLOW_ROOT}/scr/mnp-flex_get_report.sh"
EXTRACT_SCORE_SCRIPT="${WORKFLOW_ROOT}/scr/extract_report_score.R"

# Optional: check if files exist (good practice)
for file in "$MNPFLEX_BED" "$MNPFLEX_REPORT_SCRIPT" "$EXTRACT_SCORE_SCRIPT"; do
    if [ ! -e "$file" ]; then
        echo "Error: Missing file $file" >&2
        exit 1
    fi
done

# Output paths
merged_bam="/path/to/output/results/merged/${id}.merged.bam"
mods_out="/path/to/output/results/mods"
mnpflex_out="/path/to/output/results/mnpflex"

# Find all BAMs in the folder
bam_files=($(find "$bam_dir" -name "*.bam"))
if [[ ${#bam_files[@]} -eq 0 ]]; then
    echo "No BAM files found!"
    exit 1
fi
echo "Found ${#bam_files[@]} BAM files for ID $id"

# Merge BAMs
echo "Merging BAM files..."
mkdir -p "$(dirname "$merged_bam")"
samtools merge -@ "$modkit_threads" "$merged_bam" "${bam_files[@]}"

# Index merged BAM
echo "Indexing merged BAM..."
samtools index -@ "$modkit_threads" "$merged_bam"

# Run modkit pileup
echo "Running modkit pileup..."
mkdir -p "$mods_out"
modkit pileup \
    "$merged_bam" \
    "$mods_out/${id}.mods.bedmethyl" \
    --ref "$ref" \
    --threads "$modkit_threads"
echo "Done: ${mods_out}/${id}.mods.bedmethyl"

# Run MNPFlex preprocessing
mkdir -p "$mnpflex_out"
echo "Starting MNPFlex preprocessing..."
bash "$mnpflex_script" \
    "$mods_out/${id}.mods.bedmethyl" \
    "$MNPFLEX_BED" \
    "$mnpflex_out" \
    "${id}"
echo "Finished MNPFlex preprocessing."

# Upload to MNPFlex platform and get report
mkdir -p "$mnpflex_out"
echo "Uploading bed file and generating MNPFlex report..."
bash "$MNPFLEX_REPORT_SCRIPT" \
    "${mnpflex_out}/${id}.MNPFlex.subset.bed" \
    "${mnpflex_out}/${id}.MNPFlex.subset.pdf"
echo "Finished generating MNPFlex report."

# Extract classification scores and MGMT status from the pdf report
echo "Extracting scores from the MNPFlex report..."
Rscript "$EXTRACT_SCORE_SCRIPT" "${mnpflex_out}/${id}.MNPFlex.subset.pdf"
echo "Finished extracting scores from the report."


