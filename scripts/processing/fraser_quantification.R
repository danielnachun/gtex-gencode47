# fraser_quantification.R
library(FRASER)

# Function for FRASER quantification on a bam
perform_fraser_quantification <- function(bam_file, output_dir, output_name) {
    fds <- FraserDataSet()
    bamFile(fds) <- bam_file
    fds <- countRNAData(fds)
    saveFraserDataSet(fds, dir = output_dir, name = output_name)
}

# Main script to call the function with command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign command-line arguments to variables
bam_file_path <- args[1]
output_dir <- args[2]
output_name <- args[3]

# Call the function
perform_fraser_quantification(bam_file_path, output_dir, output_name)