library(SummarizedExperiment)
library(mutscan)
library(optparse)

# Define command line options
# Define command line options
option_list <- list(
  make_option(c("-f", "--infile_f"), type = "character", default = NULL, help = "Path to the input forward FASTQ file", metavar = "character"),
  make_option(c("-r", "--infile_r"), type = "character", default = NULL, help = "Path to the input reverse FASTQ file", metavar = "character"),
  make_option(c("-o", "--outname"), type = "character", default = "mutscan_output", help = "Base name for output files", metavar = "character")
)

# Parse command line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input
if (is.null(opt$infile_f) || is.null(opt$infile_r)) {
  stop("Both input files must be specified with the --infile_f and --infile_r options.", call. = FALSE)
}

# Use the current directory for output
output_dir <- paste0(getwd(), "/", format(Sys.Date(), "%Y-%m-%d"), "_", format(Sys.time(), "%H-%M-%S"), "/")

# Check if output directory exists and create it if it does not
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Modify output file name to include the --outname argument and date-time information
output_name <- paste0(output_dir, opt$outname, "_", format(Sys.Date(), "%Y-%m-%d"), "_", format(Sys.time(), "%H-%M-%S"))
output_file <- paste0(output_name, ".csv")

# Your existing code for processing goes here, make sure to replace the input file path with `opt$infile`
# For example, replace your `digestFastqs` call with the new input file path:
output <- digestFastqs(
  fastqForward = opt$infile_f,
  fastqReverse = opt$infile_r,
  elementsForward = "SPVS",
  elementsReverse = "SPVS",
  primerForward = "GCTGGTGAGGTTGCGGATAACGC",
  primerReverse = "TTAGGGAGTAGGGTAGTG",
  elementLengthsForward = c(-1, 23, 18, -1),
  elementLengthsReverse = c(-1, 18, 18, -1),
  maxReadLength = 6000  # Increase the maximum read length as needed
  # Add other arguments as needed
)

# The rest of your script remains the same, just ensure you use `output_file` for writing outputs

# Save the summary table as a CSV file
write.csv(output$summaryTable, file = output_file, row.names = FALSE)

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
png(filename = paste0(output_name, "_filtering_summary.png"), width = 800, height = 600)
plotFiltering(se, valueType = "reads", onlyActiveFilters = TRUE)
dev.off()

# Plot the error statistics and save as PNG
#png(filename = paste0("/home/ubuntu/AWTallyNN/tallynn/R/", format(date, "%Y-%m-%d"), "_error_stats.png"), width = 800, height = 600)
#plotErrorStats(se, valueType = "reads")
#dev.off()

# Plot the distribution of variant counts and save as PNG
png(filename = paste0(output_name, "_variant_distribution.png"), width = 800, height = 600)
plotDistributions(se, plotType = "density", pseudocount = 1)
dev.off()

# Plot the total counts per variant and save as PNG
png(filename = paste0(output_name, "_variant_totals.png"), width = 800, height = 600)
plotTotals(se)
dev.off()

