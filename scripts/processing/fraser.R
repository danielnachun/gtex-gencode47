#!/usr/bin/env Rscript

library(argparser)
library(S4Vectors)
library(FRASER)
library(BSgenome)
library(readr)
library(stringr)
library(BiocManager)

# parse arguments
parser <- argparser::arg_parser("Script to count split reads with FRASER") |>
    argparser::add_argument("--bam_file", help = "Path to BAM file") |>
    argparser::add_argument("--sample_id", help = "Name of sample") |>
    argparser::add_argument("--output_directory", help = "Path to output directory") |>
    argparser::add_argument("--working_directory", help = "Path to working directory")

parsed_args <- argparser::parse_args(parser)
bam_file <- parsed_args$bam_file
sample_id <- parsed_args$sample_id
output_directory <- parsed_args$output_directory
working_directory <- parsed_args$working_directory

# get the genome
install("BSgenome.Hsapiens.UCSC.hg38")
genome <- BSgenome::getBSgenome("hg38")

# load in the sample bam
sample_annot <- S4Vectors::DataFrame(sampleID = sample_id, bamFile = bam_file)
fraser_object <- FRASER::FraserDataSet(colData = sample_annot, workingDir = working_directory, name = stringr::str_c("raw-local-", sample_id))
FRASER::pairedEnd(fraser_object) <- TRUE
FRASER::strandSpecific(fraser_object) <- FALSE

# count the reads
split_reads <- FRASER::countSplitReads(sampleID = sample_id, fds = fraser_object, NcpuPerSample = 1, recount = TRUE, keepNonStandardChromosomes = FALSE, genome = genome)

# write out the rds
readr::write_rds(split_reads, stringr::str_c(output_directory, "/", sample_id, ".split_reads.rds"))
