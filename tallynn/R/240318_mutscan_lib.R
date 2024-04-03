library(SummarizedExperiment)
library(mutscan)

output <- digestFastqs(
  fastqForward = "/home/ubuntu/AWTallyNN/tallynn/python/240328_Cloned_library",
  elementsForward = "SPVS",
  primerForward = "GCTGGTGAGGTTGCGGATAACGC",
  elementLengthsForward = c(-1, 23, 18, -1),
  maxReadLength = 5000  # Increase the maximum read length as needed
  # Add other arguments as needed
)

date <- Sys.Date()

# Specify the output file path and name for the CSV file, including the date
output_dir <- "/home/ubuntu/AWTallyNN2/AWTallyNN/tallynn/R/" + format(date, "%Y-%m-%d") + "/"
output_file <- paste0(output_dir, "mutscan_output_", format(date, "%Y-%m-%d"), ".csv")

# Extract the summary table from the output object
summary_table <- output$summaryTable

# Save the summary table as a CSV file
write.csv(summary_table, file = output_file, row.names = FALSE)

# Generate SummarizedExperiment object
se <- summarizeExperiment(
  x = list(sample = output),
  coldata = data.frame(Name = "sample", Condition = "condition"),
  countType = "reads"  # Set countType to "reads" instead of "umis"
)

# Access count matrix, sample annotations, and variant annotations
head(assay(se, "counts"))
Matrix::colSums(assay(se, "counts"))
head(rowData(se))
colData(se)

# Check count type (reads or UMIs)
metadata(se)$countType

# Create plots using mutscan functions and save them as PNG files

# Plot the filtering summary and save as PNG
png(filename = paste0(output_dir, format(date, "%Y-%m-%d"), "_filtering_summary.png"), width = 800, height = 600)
plotFiltering(se, valueType = "reads", onlyActiveFilters = TRUE)
dev.off()

# Plot the error statistics and save as PNG
#png(filename = paste0("/home/ubuntu/AWTallyNN/tallynn/R/", format(date, "%Y-%m-%d"), "_error_stats.png"), width = 800, height = 600)
#plotErrorStats(se, valueType = "reads")
#dev.off()

# Plot the distribution of variant counts and save as PNG
png(filename = paste0(output_dir, format(date, "%Y-%m-%d"), "_variant_distribution.png"), width = 800, height = 600)
plotDistributions(se, plotType = "density", pseudocount = 1)
dev.off()

# Plot the total counts per variant and save as PNG
png(filename = paste0(output_dir, format(date, "%Y-%m-%d"), "_variant_totals.png"), width = 800, height = 600)
plotTotals(se)
dev.off()

