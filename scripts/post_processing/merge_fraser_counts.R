library(FRASER)
library(furrr)
library(GenomicRanges)
library(tidyverse)
library(magrittr)

rds_files <- list.files("../../output/all_tissues_quantifications/fraser/", full.names = TRUE)
gtex_tissues <- readr::read_tsv("../../data/other_references/v10/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
output_directory <- "../../output/fraser_aggregated"

dir.create(output_directory, showWarnings = FALSE)

file_names_cleaned <- stringr::str_remove_all(rds_files, "\\.split_reads.rds") |> basename()
gtex_tissues_filter <- dplyr::filter(gtex_tissues, magrittr::is_in(SAMPID, file_names_cleaned))
gtex_tissues_filter$tissue_simplified  <- stringr::str_remove_all(gtex_tissues_filter$SMTSD, "Cells - ") |> stringr::str_remove_all(" - .*$")

aggregate_spliced_reads <- function(grouped_rows, grouping_df, rds_files, output_directory) {
    tissue_rds <- stringr::str_subset(rds_files, stringr::str_c(grouped_rows$SAMPID, collapse = "|"))
    sample_ids <- stringr::str_remove_all(tissue_rds, "\\.split_reads.rds") |> basename()

    future::plan(multisession, workers = 16)
    tissue_list_collapse <- furrr::future_map(tissue_rds, read_rds, .progress = TRUE) |>
        GenomicRanges::GRangesList() |> magrittr::set_names(sample_ids) |>
        unlist() |> unique() |> sort()
    tissue_collapse_name <- stringr::str_c(output_directory, "/", grouping_df$tissue_simplified, ".rds")
    readr::write_rds(tissue_list_collapse, tissue_collapse_name)
    grouping_df$tissue_simplified
}

tissue_list_collapse <- dplyr::group_by(gtex_tissues_filter, tissue_simplified) |>
    dplyr::group_map(aggregate_spliced_reads, rds_files, output_directory)

# Commented out because this is too memory intensive right now
# get_overlaps <- function(sample_name, gr, ranges, output_directory) {
#     sample_count <- integer(length(ranges))
#
#     overlaps <- findOverlaps(gr, ranges, type = "equal")
#     sample_count[overlaps@to] <- mcols(gr)$count
#     file_path <- str_c(output_directory, "/", sample_name, ".tsv.gz")
#     sample_df <- as_tibble(ranges) |> select(seqnames, start, end, strand)
#     sample_df$counts <- sample_count
#     write_tsv(sample_df, file_path)
# }
#
# sample_counts <- future_map2(names(tissue_list_gr), tissue_list_gr, get_overlaps, tissue_list_collapse, output_directory, .progress = TRUE)

temp_tissue <- filter(gtex_tissues_filter, SMTS == "Brain") |> slice(1)
sample_id <- temp_tissue$SAMPID
bam_file <- str_c("/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/all_tissues_quantifications/genome_bam/", sample_id, ".Aligned.sortedByCoord.out.patched.v11md.bam")
working_directory <- "/tmp/test"

# load in the sample bam
sample_annot <- S4Vectors::DataFrame(sampleID = sample_id, bamFile = bam_file)
fraser_object <- FRASER::FraserDataSet(colData = sample_annot, workingDir = working_directory, name = stringr::str_c("raw-local-", sample_id))
FRASER::pairedEnd(fraser_object) <- TRUE
FRASER::strandSpecific(fraser_object) <- FALSE
test_ranges <- read_rds("/oak/stanford/groups/smontgom/dnachun/data/gtex/v10/output/fraser_aggregated/Brain.rds")

# count the reads
split_reads <- countNonSplicedReads(sampleID = sample_id, fds = fraser_object, NcpuPerSample = 1, recount = TRUE, spliceSiteCoords = test_ranges)
