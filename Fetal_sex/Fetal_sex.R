#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

# ============================================================
# NEWPAT - Fetal sex determination
# ============================================================

option_list <- list(

  # Required inputs

  make_option("--autosomes", type = "character",
              help = "Autosomal DNPcall output from maternal cfDNA"),

  make_option("--X", type = "character",
              help = "X chromosome DNPcall output from maternal cfDNA"),

  make_option("--Y", type = "character",
              help = "Y chromosome SNPcall output from maternal cfDNA"),

  make_option("--out", type = "character",
              help = "Output directory"),

  # Marker and fetal read selection

  make_option("--max_other_fraction",
              type = "double",
              default = 0.10,
              help = "Maximum Other / total_count fraction [default: %default]"),

  make_option("--maternal_hom_threshold",
              type = "double",
              default = 0.90,
              help = "Maternal homozygosity threshold [default: %default]"),

  make_option("--include_zero_coverage_Y",
              type = "character",
              default = "TRUE",
              help = "Include zero REF+ALT Y markers: TRUE or FALSE [default: %default]"),


  make_option("--min_fetal_reads", type = "integer",
              default = NULL,
              help = "Minimum non-maternal reads"),

  make_option("--min_fetal_percent", type = "double",
              default = NULL,
              help = "Minimum percentage of non-maternal reads"),

  # Fetal sex classification thresholds

  make_option("--A_threshold", type = "double",
              default = 0.5),

  make_option("--B_female_max", type = "double",
              default = 0.1),

  make_option("--B_male_min", type = "double",
              default = 1),

  make_option("--C_female_max", type = "double",
              default = 0.1),

  make_option("--C_male_min", type = "double",
              default = 2)
)

opt <- parse_args(
  OptionParser(option_list = option_list)
)

# ============================================================
# 1. INPUT VALIDATION
# ============================================================

required <- c("autosomes", "X", "Y", "out")

for (arg in required) {
  if (is.null(opt[[arg]]) || !nzchar(opt[[arg]])) {
    stop("Missing required argument: --", arg)
  }
}

if (is.null(opt$min_fetal_reads) &&
    is.null(opt$min_fetal_percent)) {

  opt$min_fetal_reads <- 2L
}

if (!is.null(opt$min_fetal_reads)) {
  if (is.na(opt$min_fetal_reads) ||
      opt$min_fetal_reads < 1) {
    stop("--min_fetal_reads must be >= 1")
  }
}

if (!is.null(opt$min_fetal_percent)) {
  if (!is.finite(opt$min_fetal_percent) ||
      opt$min_fetal_percent <= 0 ||
      opt$min_fetal_percent > 100) {
    stop("--min_fetal_percent must be > 0 and <= 100")
  }
}

if (!is.finite(opt$max_other_fraction) ||
    opt$max_other_fraction < 0 ||
    opt$max_other_fraction > 1) {
  stop("--max_other_fraction must be between 0 and 1")
}

if (!is.finite(opt$maternal_hom_threshold) ||
    opt$maternal_hom_threshold <= 0.5 ||
    opt$maternal_hom_threshold > 1) {
  stop("--maternal_hom_threshold must be > 0.5 and <= 1")
}

opt$include_zero_coverage_Y <-
  toupper(opt$include_zero_coverage_Y)

if (!opt$include_zero_coverage_Y %in% c("TRUE", "FALSE")) {
  stop("--include_zero_coverage_Y must be TRUE or FALSE")
}

opt$include_zero_coverage_Y <-
  opt$include_zero_coverage_Y == "TRUE"

thresholds <- c(
  opt$A_threshold,
  opt$B_female_max,
  opt$B_male_min,
  opt$C_female_max,
  opt$C_male_min
)

if (any(!is.finite(thresholds)) ||
    any(thresholds < 0) ||
    opt$A_threshold <= 0 ||
    opt$B_female_max >= opt$B_male_min ||
    opt$C_female_max >= opt$C_male_min) {
  stop("Invalid fetal sex classification thresholds")
}

if (file.exists(opt$out)) {
  stop("Output path already exists: ", opt$out)
}

# ============================================================
# 2. FUNCTIONS FOR READING INPUT FILES
# ============================================================

marker_key <- c("chr", "pos1", "pos2", "REF", "ALT")

read_table <- function(path, required_columns) {

  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  x <- read.delim(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  missing_columns <- setdiff(
    required_columns,
    names(x)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in ", path, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  return(x)
}

read_calling <- function(path, marker_type, region) {

  if (marker_type == "DNP") {
    key <- marker_key
  } else {
    key <- c("chr", "pos", "REF", "ALT")
  }

  count_columns <- c(
    "N_REF", "N_ALT", "N_Other", "total_count"
  )

  x <- read_table(
    path,
    c(key, count_columns)
  )

  if (nrow(x) == 0) {
    stop("Empty calling file: ", path)
  }

  if (anyNA(x[, key, drop = FALSE])) {
    stop("Missing marker identifiers in: ", path)
  }

  if (anyDuplicated(x[, key, drop = FALSE])) {
    stop("Duplicated markers in: ", path)
  }

  for (column in count_columns) {

    values <- x[[column]]

    if (!is.numeric(values) ||
        any(!is.finite(values)) ||
        any(values < 0) ||
        any(values != floor(values))) {

      stop("Invalid read counts in: ", path)
    }
  }

  if (any(x$total_count <
          x$N_REF + x$N_ALT + x$N_Other)) {
    stop("Inconsistent total_count in: ", path)
  }

  chr <- sub(
    "^chr", "", as.character(x$chr),
    ignore.case = TRUE
  )

  if (region == "autosomes" &&
      any(!chr %in% as.character(1:22))) {
    stop("Non-autosomal markers in: ", path)
  }

  if (region == "X" && any(chr != "X")) {
    stop("Non-X markers in: ", path)
  }

  if (region == "Y" && any(chr != "Y")) {
    stop("Non-Y markers in: ", path)
  }

  return(x)
}

write_output <- function(x, filename) {

  write.table(
    x,
    file.path(opt$out, filename),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
}

# ============================================================
# 3. READ MATERNAL cfDNA CALLING OUTPUTS
# ============================================================

autosomes <- read_calling(
  opt$autosomes, "DNP", "autosomes"
)

X <- read_calling(
  opt$X, "DNP", "X"
)

Y <- read_calling(
  opt$Y, "SNP", "Y"
)

# X marker specificity must be validated upstream.

# ============================================================
# 4. MATERNAL INFORMATIVE SITE SELECTION
# ============================================================

compute_madre_info <- function(x) {

  info <- x %>%
    mutate(
      ratio_Ref_mother =
        N_REF / (N_REF + N_ALT),

      ratio_Other =
        N_Other / total_count
    ) %>%
    filter(
      N_REF + N_ALT > 0,
      total_count > 0,
      ratio_Other <= opt$max_other_fraction,
      ratio_Ref_mother >= opt$maternal_hom_threshold |
      ratio_Ref_mother <= round(1 - opt$maternal_hom_threshold, 15)
    ) %>%
    mutate(

      fetal_reads = ifelse(
        ratio_Ref_mother >= opt$maternal_hom_threshold,
        N_ALT,
        N_REF
      ),

      fetal_percent =
        100 * fetal_reads /
        (N_REF + N_ALT)
    )

  if (!is.null(opt$min_fetal_reads)) {

    info <- info %>%
      filter(
        fetal_reads >= opt$min_fetal_reads
      )
  }

  if (!is.null(opt$min_fetal_percent)) {

    info <- info %>%
      filter(
        fetal_percent >= opt$min_fetal_percent
      )
  }

  return(info)
}

info_A <- compute_madre_info(autosomes)

info_X <- compute_madre_info(X)

# ============================================================
# 5. Y CHROMOSOME READS
# ============================================================

info_Y <- Y %>%
  mutate(
    Y_reads = N_REF + N_ALT,
    Covered = Y_reads > 0,
    Included = opt$include_zero_coverage_Y | Covered
  )

# ============================================================
# 6. CALCULATE MFF
# ============================================================

MFF_A <- if (nrow(info_A) > 0) {

  mean(info_A$fetal_reads)

} else {

  NA_real_
}

MFF_X <- if (nrow(info_X) > 0) {

  mean(info_X$fetal_reads)

} else {

  NA_real_
}

MFF_Y <- if (any(info_Y$Included)) {
  mean(info_Y$Y_reads[info_Y$Included])
} else {
  0
}

# ============================================================
# 7. CALCULATE A, B AND C RATIOS
# ============================================================

safe_ratio <- function(numerator, denominator) {

  if (is.na(numerator) ||
      is.na(denominator)) {

    return(NA_real_)
  }

  if (denominator == 0) {

    if (numerator > 0) {
      return(Inf)
    }

    return(NA_real_)
  }

  return(numerator / denominator)
}

A <- safe_ratio(MFF_X, MFF_A)

B <- safe_ratio(MFF_Y, MFF_A)

C <- safe_ratio(MFF_Y, MFF_X)

# ============================================================
# 8. FETAL SEX CLASSIFICATION
# ============================================================

classify_sex <- function(A, B, C) {

  if (any(is.na(c(A, B, C)))) {

    return("Undetermined")
  }

  if (A >= opt$A_threshold &&
      B <= opt$B_female_max &&
      C <= opt$C_female_max) {

    return("F")
  }

  if (A < opt$A_threshold &&
      B >= opt$B_male_min &&
      C >= opt$C_male_min) {

    return("M")
  }

  return("Undetermined")
}

fetal_sex <- classify_sex(A, B, C)

# ============================================================
# 9. SUMMARY
# ============================================================

summary_result <- data.frame(

  Sample = sub(
    "^full_output_", "",
    tools::file_path_sans_ext(
      basename(opt$autosomes)
    )
  ),

  N_A_original = nrow(autosomes),

  N_A_sites = nrow(info_A),

  N_X_original = nrow(X),

  N_X_sites = nrow(info_X),

  N_Y_sites = nrow(info_Y),

  N_Y_covered_sites = sum(info_Y$Covered),

  N_Y_sites_used = sum(info_Y$Included),

  Fetal_reads_A = sum(info_A$fetal_reads),

  Fetal_reads_X = sum(info_X$fetal_reads),

  Reads_Y = sum(info_Y$Y_reads),

  MFF_A = MFF_A,

  MFF_X = MFF_X,

  MFF_Y = MFF_Y,

  Ratio_A = A,

  Ratio_B = B,

  Ratio_C = C,

  Fetal_sex = fetal_sex,

  Max_other_fraction = opt$max_other_fraction,

  Maternal_hom_threshold = opt$maternal_hom_threshold,

  Include_zero_coverage_Y = opt$include_zero_coverage_Y,

  Min_fetal_reads = if (
    is.null(opt$min_fetal_reads)
  ) NA_integer_ else opt$min_fetal_reads,

  Min_fetal_percent = if (
    is.null(opt$min_fetal_percent)
  ) NA_real_ else opt$min_fetal_percent,

  A_threshold = opt$A_threshold,

  B_female_max = opt$B_female_max,

  B_male_min = opt$B_male_min,

  C_female_max = opt$C_female_max,

  C_male_min = opt$C_male_min
)

# ============================================================
# 10. SAVE OUTPUT FILES
# ============================================================

dir.create(
  opt$out,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(opt$out)) {
  stop("Cannot create output directory")
}

write_output(
  summary_result,
  "Fetal_sex_summary.txt"
)

write_output(
  info_A,
  "Autosomal_informative_sites.txt"
)

write_output(
  info_X,
  "X_informative_sites.txt"
)

write_output(
  info_Y,
  "Y_sites.txt"
)

# ============================================================
# 11. TERMINAL OUTPUT
# ============================================================

cat("\n====================================\n")
cat("NEWPAT - FETAL SEX DETERMINATION\n")
cat("====================================\n\n")

cat("Sample:", summary_result$Sample, "\n")

cat("Autosomal informative sites:",
    nrow(info_A), "\n")

cat("X informative sites:",
    nrow(info_X), "\n")

cat("Y covered sites:",
    sum(info_Y$Covered), "/", nrow(info_Y), "\n\n")

cat("Y sites used:", sum(info_Y$Included), "\n\n")

cat("MFF_A:", MFF_A, "\n")
cat("MFF_X:", MFF_X, "\n")
cat("MFF_Y:", MFF_Y, "\n\n")

cat("Ratio A (X/autosomes):", A, "\n")
cat("Ratio B (Y/autosomes):", B, "\n")
cat("Ratio C (Y/X):", C, "\n\n")

cat("Fetal sex:", fetal_sex, "\n")

cat("\nResults saved in:", opt$out, "\n")

