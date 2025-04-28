library(pdftools)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("No input file provided. Usage: Rscript script.R <report_file.pdf>")
}

report_file <- args[1]
output_file <- sub("pdf$", "txt", report_file)

results <- list()

# Read PDF
pdf_text_data <- pdf_text(report_file)

# --- Extract Methylation classification scores ---

i <- grep("^Methylation classification\\s+Methylation classification", pdf_text_data)
text_page <- pdf_text_data[i]

# Split into lines
lines <- strsplit(text_page, "\n")[[1]]

# Initialize empty data.frame
class_scores <- NULL

for (j in c("Methylation superfamily", "Methylation family", "Methylation class", "Methylation subclass")) {
  score_line <- grep(paste0(".+", j, ".+[0-9\\.]+$"), lines, value = TRUE)
  
  # Clean and split
  score_parts <- strsplit(score_line, "\\s{2,}")[[1]] # split by 2+ spaces
  score_parts <- score_parts[score_parts != ""]        # remove empty parts
  
  class_scores <- rbind(class_scores, data.frame(
    name = j,
    class = score_parts[2],
    score = as.numeric(score_parts[3]),
    stringsAsFactors = FALSE
  ))
}

results[["class_scores"]] <- class_scores

# --- Extract MGMT promoter methylation info ---

i <- grep("Methylation classification\\s+MGMT promoter methylation", pdf_text_data)
text_page <- pdf_text_data[i]

# Extract fields
num_sites <- str_match(text_page, "Number of sites used =\\s*(\\d+)")[,2]
avg_meth <- str_match(text_page, "Average methylation =\\s*(\\d+)%")[,2]
mgmt_status <- str_match(text_page, "Predicted MGMT promoter status=\\s*(.*)")[,2]

# Store
results[["MGMT"]] <- c(
  "Number of sites used" = as.integer(num_sites),
  "Average methylation" = as.numeric(avg_meth),
  "Predicted MGMT promoter status" = mgmt_status
)

# --- Save results ---
dput(results, file = output_file)

# Later you can read by:
# restored_results <- dget(output_file)
