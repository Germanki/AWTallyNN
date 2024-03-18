library(mutscan)

output <- digestFastqs(
  fastqForward = "/home/ubuntu/AWTallyNN/tallynn/python/240318_IDTGblock_Flipped_2",
  elementsForward = "SPVS",
  primerForward = "GCTGGTGAGGTTGCGGATAACG",
  elementLengthsForward = c(-1, 22, 18, -1),
  maxReadLength = 2000  # Increase the maximum read length as needed
  # Add other arguments as needed
)

date <- Sys.Date()

# Specify the output file path and name for the CSV file, including the date
output_file <- paste0("/home/ubuntu/AWTallyNN/tallynn/R/mutscan_output_", format(date, "%Y-%m-%d"), ".csv")

# Extract the summary table from the output object
summary_table <- output$summaryTable

# Save the summary table as a CSV file
write.csv(summary_table, file = output_file, row.names = FALSE)

# Create plots using mutscan functions and save them as PNG files

# Plot the filtering summary and save as PNG
png(filename = "/home/ubuntu/AWTallyNN/tallynn/R/", format(date, "%Y-%m-%d"), "_filtering_summary.png", width = 800, height = 600)
plotFiltering(output, valueType = "reads", onlyActiveFilters = TRUE)
dev.off()

# Plot the error statistics and save as PNG
png(filename = "/home/ubuntu/AWTallyNN/tallynn/R/", format(date, "%Y-%m-%d"), "_error_stats.png", width = 800, height = 600)
plotErrorStats(output)
dev.off()

# Plot the distribution of variant counts and save as PNG
png(filename = "/home/ubuntu/AWTallyNN/tallynn/R/", format(date, "%Y-%m-%d"), "_variant_distribution.png", width = 800, height = 600)
plotDistributions(output, plotType = "density", pseudocount = 1)
dev.off()

# Plot the total counts per variant and save as PNG
png(filename = "/home/ubuntu/AWTallyNN/tallynn/R/", format(date, "%Y-%m-%d"), "_variant_totals.png", width = 800, height = 600)
plotTotals(output)
dev.off()

