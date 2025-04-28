#!/bin/bash

set -euo pipefail

# Inputs
bam_dir="/path/to/folder/of/aligned/bams"             # Directory containing your BAM files
ref="/path/to/GRCh38.p14.genome.fa"
id="sample_name"
modkit_threads=8
mnpflex_script="/path/to/Rapid-CNS2_nf/scr/mnp-flex_preprocessing.sh"
mnpflex_bed="/path/to/Rapid-CNS2_nf/data/MNP-flex.bed"
mnpflex_report_script="/path/to/Rapid-CNS2_nf/scr/mnp-flex_get_report.sh"
extract_score_script="/path/to/Rapid-CNS2_nf/scr/extract_report_score.R"

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
    "$mnpflex_bed" \
    "$mnpflex_out" \
    "${id}"
echo "Finished MNPFlex preprocessing."

# Upload to MNPFlex platform and get report
mkdir -p "$mnpflex_out"
echo "Uploading bed file and generating MNPFlex report..."
bash "$mnpflex_report_script" \
    "${mnpflex_out}/${id}.MNPFlex.subset.bed" \
    "${mnpflex_out}/${id}.MNPFlex.subset.pdf"
echo "Finished generating MNPFlex report."

# Extract classification scores and MGMT status from the pdf report
echo "Extracting scores from the MNPFlex report..."
Rscript "$extract_score_script" "${mnpflex_out}/${id}.MNPFlex.subset.pdf"
echo "Finished extracting scores from the report."


