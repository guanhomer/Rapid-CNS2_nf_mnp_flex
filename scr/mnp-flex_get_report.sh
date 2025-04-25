#!/bin/bash

# --- File settings ---
# input_sample_path="/Users/rkwan/Downloads/mnp_flex_sample_data.bed"
# output_report_path="/Users/rkwan/Downloads/mnp_flex_sample_data.pdf"
input_sample_path=$1
output_report_path=$2

# --- Assign a random sample name ---
sample_name="sample_$(date +%Y%m%d%H%M%S)"

# --- Get token ---
echo "Getting authentication token..."
tocken=$(curl -s -X POST https://mnp-flex.org/api/v1/auth/token \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "username=demo_user" \
  -d "password=KLrHgx7HMC" \
  -d "grant_type=" \
  -d "scope=" \
  -d "client_id=" \
  -d "client_secret=" | jq -r '.access_token')

echo "Token received."

# --- Upload sample ---
echo "Uploading sample: $sample_name ..."
upload_response=$(curl -s -X PUT "https://mnp-flex.org/api/v1/samples/?sample_name=${sample_name}&disclaimer_confirmed=true" \
  -H "Authorization: Bearer ${tocken}" \
  -F "file=@${input_sample_path}")

echo "Upload response: $upload_response"

# --- Extract status and sample ID from upload response ---
analysis_status=$(echo "$upload_response" | jq -r '.bed_file_sample.analysis_status')
sample_id=$(echo "$upload_response" | jq -r '.id')

echo "Extracted: Sample ID = $sample_id, Analysis Status = $analysis_status"

# --- Check analysis status ---
if [ "$analysis_status" == "done" ]; then
  echo "Analysis complete immediately. Downloading report..."

  curl -s -X GET "https://mnp-flex.org/api/v1/samples/download_sample_result/${sample_id}" \
    -H "Authorization: Bearer ${tocken}" \
    --output "${output_report_path}"

  echo "Report downloaded to ${output_report_path}"
  exit 0
else
  echo "Analysis not done yet. Start polling..."
  # (Optional: fallback to polling logic here if needed)
fi