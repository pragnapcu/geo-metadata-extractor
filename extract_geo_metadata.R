#!/usr/bin/env Rscript
# GEOquery sample metadata extraction for local GEO SOFT family files.
# Usage: Rscript extract_geo_metadata.R GSE178994_family.soft.gz output_dir

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: Rscript extract_geo_metadata.R <family.soft[.gz]> [output_dir]", call. = FALSE)
}
input_file <- args[[1L]]
out_dir <- if (length(args) == 2L) args[[2L]] else "geo_metadata_output"
if (!file.exists(input_file)) stop("Input file not found: ", input_file, call. = FALSE)
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop("GEOquery is not installed. In R run: BiocManager::install('GEOquery')", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading local GEO SOFT file: ", input_file)
gse <- GEOquery::getGEO(filename = input_file, GSEMatrix = FALSE)
samples <- GEOquery::GSMList(gse)
if (length(samples) == 0L) stop("No GSM samples were found.", call. = FALSE)
series_meta <- GEOquery::Meta(gse)
gse_id <- as.character(series_meta$geo_accession)[1L]
if (is.na(gse_id) || !nzchar(gse_id)) {
  gse_id <- sub(".*(GSE[0-9]+).*", "\\1", basename(input_file))
}
if (!grepl("^GSE[0-9]+$", gse_id)) gse_id <- "UNKNOWN_GSE"

# Keep missing and empty metadata distinct from observations.
get_field <- function(m, key) {
  x <- m[[key]]
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0L) return(NA_character_)
  paste(x, collapse = " | ")
}

# Extract every value with a given characteristic key; never silently take only the first.
get_characteristic <- function(m, key) {
  x <- m[["characteristics_ch1"]]
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  x <- as.character(x)
  pattern <- paste0("^\\s*", key, "\\s*:")
  selected <- grep(pattern, x, ignore.case = TRUE, value = TRUE)
  if (length(selected) == 0L) return(NA_character_)
  values <- trimws(sub("^[^:]*:", "", selected))
  values <- values[nzchar(values)]
  if (length(values) == 0L) return(NA_character_)
  paste(unique(values), collapse = " | ")
}

raw_list <- lapply(names(samples), function(gsm) {
  m <- GEOquery::Meta(samples[[gsm]])
  if (length(m) == 0L) return(NULL)
  # Rows per ORIGINAL metadata value, including repeated characteristics_ch1.
  data.frame(
    gse = rep(gse_id, sum(lengths(m))),
    gsm = rep(gsm, sum(lengths(m))),
    field = rep(names(m), lengths(m)),
    value = as.character(unlist(m, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
})
raw_df <- do.call(rbind, raw_list)

all_fields <- sort(unique(raw_df$field))
field_counts <- data.frame(
  field = all_fields,
  sample_count = vapply(all_fields, function(field) {
    sum(vapply(samples, function(s) field %in% names(GEOquery::Meta(s)), logical(1)))
  }, integer(1)),
  stringsAsFactors = FALSE
)

structured <- do.call(rbind, lapply(names(samples), function(gsm) {
  m <- GEOquery::Meta(samples[[gsm]])
  data.frame(
    gse = gse_id,
    gsm = gsm,
    title = get_field(m, "title"),
    description = get_field(m, "description"),
    disease = get_characteristic(m, "disease state"),
    cell_type = get_characteristic(m, "cell type"),
    stimulation = get_characteristic(m, "stimulation"),
    tissue = get_field(m, "source_name_ch1"),
    organism = get_field(m, "organism_ch1"),
    library_strategy = get_field(m, "library_strategy"),
    treatment_protocol = get_field(m, "treatment_protocol_ch1"),
    relation = get_field(m, "relation"),
    characteristics_raw = get_field(m, "characteristics_ch1"),
    stringsAsFactors = FALSE
  )
}))

write.csv(raw_df, file.path(out_dir, paste0(gse_id, "_raw_metadata.csv")), row.names = FALSE, na = "")
write.csv(field_counts, file.path(out_dir, paste0(gse_id, "_metadata_field_counts.csv")), row.names = FALSE, na = "")
write.csv(structured, file.path(out_dir, paste0(gse_id, "_structured_metadata.csv")), row.names = FALSE, na = "")

# Data summaries include NA as an explicit category, not an inferred clinical label.
summary_fields <- c("disease", "cell_type", "stimulation", "tissue")
summary_rows <- lapply(summary_fields, function(field) {
  tab <- table(structured[[field]], useNA = "ifany")
  data.frame(attribute = field, value = names(tab), count = as.integer(tab), stringsAsFactors = FALSE)
})
summary_df <- do.call(rbind, summary_rows)
summary_df$value[is.na(summary_df$value)] <- "<MISSING>"
write.csv(summary_df, file.path(out_dir, paste0(gse_id, "_sample_summary.csv")), row.names = FALSE)

message("GSE: ", gse_id)
message("GSM samples: ", length(samples))
message("Platform records: ", length(GEOquery::GPLList(gse)))
message("Raw metadata rows: ", nrow(raw_df))
message("Unique metadata fields: ", length(all_fields))
message("Output directory: ", normalizePath(out_dir))
print(summary_df, row.names = FALSE)
message("Finished. No LLM harmonization is performed by this script.")
