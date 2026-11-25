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

  # Optional inputs

  make_option("--male_controls", type = "character",
              default = NULL,
              help = "Text file containing paths to male X DNPcall outputs"),

  make_option("--exclude_X", type = "character",
              default = NULL,
              help = "Table of X markers to exclude"),

  make_option("--exclude_untested_X",
              action = "store_true",
              default = FALSE,
              help = "Exclude X markers not evaluated in male controls"),

  # Fetal read selection

  make_option("--min_fetal_reads", type = "integer",
              default = NULL,
              help = "Minimum non-maternal reads"),

  make_option("--min_fetal_percent", type = "double",
              default = NULL,
              help = "Minimum percentage of non-maternal reads"),

  # X specificity control

  make_option("--min_male_depth", type = "integer",
              default = 100),

  make_option("--max_minor_fraction", type = "double",
              default = 0.01),

  make_option("--min_males_failed", type = "integer",
              default = 2),

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

if (is.na(opt$min_male_depth) ||
    opt$min_male_depth < 1) {
  stop("--min_male_depth must be >= 1")
}

if (!is.finite(opt$max_minor_fraction) ||
    opt$max_minor_fraction < 0 ||
    opt$max_minor_fraction >= 0.5) {
  stop("--max_minor_fraction must be >= 0 and < 0.5")
}

if (is.na(opt$min_males_failed) ||
    opt$min_males_failed < 1) {
  stop("--min_males_failed must be >= 1")
}

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

N_X_original <- nrow(X)

X_panel <- X[, marker_key, drop = FALSE]

excluded_X <- X_panel[0, ]

qc_summary <- NULL

qc_status <- "Not_performed"

# ============================================================
# 4. OPTIONAL X SPECIFICITY CONTROL
# ============================================================

if (!is.null(opt$male_controls)) {

  if (!file.exists(opt$male_controls)) {
    stop("Male control list not found")
  }

  male_files <- trimws(
    readLines(opt$male_controls, warn = FALSE)
  )

  male_files <- male_files[
    nzchar(male_files) &
    !startsWith(male_files, "#")
  ]

  if (length(male_files) == 0) {
    stop("Male control list is empty")
  }

  male_files <- normalizePath(
    male_files,
    mustWork = TRUE
  )

  if (anyDuplicated(male_files)) {
    stop("Duplicated male control files")
  }

  if (length(male_files) < opt$min_males_failed) {
    stop(
      "Insufficient male controls. ",
      "Adjust --min_males_failed or add controls."
    )
  }

  qc <- bind_rows(lapply(
    seq_along(male_files),
    function(i) {

      male_X <- read_calling(
        male_files[i], "DNP", "X"
      )

      male_X <- inner_join(
        male_X,
        X_panel,
        by = marker_key
      )

      male_X %>%
        mutate(
          Sample = i,
          Depth = N_REF + N_ALT,
          Minor_reads = pmin(N_REF, N_ALT),
          Minor_fraction = ifelse(
            Depth > 0,
            Minor_reads / Depth,
            NA_real_
          )
        ) %>%
        filter(
          Depth >= opt$min_male_depth
        ) %>%
        mutate(
          Failed = Minor_fraction >
            opt$max_minor_fraction
        )
    }
  ))

  qc_summary <- qc %>%
    group_by(
      chr, pos1, pos2, REF, ALT
    ) %>%
    summarise(
      N_males = n_distinct(Sample),
      N_failed = sum(Failed),
      Mean_minor_percent =
        100 * mean(Minor_fraction),
      Max_minor_percent =
        100 * max(Minor_fraction),
      .groups = "drop"
    )

  qc_summary <- X_panel %>%
    left_join(
      qc_summary,
      by = marker_key
    ) %>%
    mutate(
      N_males = coalesce(N_males, 0L),
      N_failed = coalesce(N_failed, 0L),
      QC_status = case_when(

        N_failed >= opt$min_males_failed ~
          "Excluded",

        N_males < opt$min_males_failed ~
          "Not_evaluated",

        TRUE ~ "Retained"
      )
    )

  excluded_X <- qc_summary %>%
    filter(QC_status == "Excluded") %>%
    select(all_of(marker_key))

  qc_status <- "Male_controls"

  if (any(qc_summary$QC_status == "Not_evaluated")) {

    warning(
      "Some X markers could not be evaluated ",
      "in enough male controls"
    )

    qc_status <- "Male_controls_partial"
  }
}

# ============================================================
# 5. OPTIONAL PREDEFINED X EXCLUSION LIST
# ============================================================

if (!is.null(opt$exclude_X)) {

  external_exclusions <- read_table(
    opt$exclude_X,
    marker_key
  )

  external_exclusions <- external_exclusions %>%
    select(all_of(marker_key))

  unmatched <- anti_join(
    external_exclusions,
    X_panel,
    by = marker_key
  )

  if (nrow(unmatched) > 0) {
    warning(
      nrow(unmatched),
      " exclusion markers do not match the X panel"
    )
  }

  external_exclusions <- semi_join(
    external_exclusions,
    X_panel,
    by = marker_key
  )

  excluded_X <- bind_rows(
    excluded_X,
    external_exclusions
  ) %>%
    distinct()

  if (is.null(opt$male_controls)) {

    qc_status <- "Predefined_exclusions"

  } else {

    qc_status <- paste0(
      qc_status,
      "_and_predefined_exclusions"
    )
  }
}

if (qc_status == "Not_performed") {

  warning(
    "X marker specificity was not evaluated. ",
    "All X markers will be retained."
  )
}

# Optional exclusion of unevaluated X markers

if (opt$exclude_untested_X) {

  if (is.null(qc_summary)) {
    stop("--exclude_untested_X requires --male_controls")
  }

  untested_X <- qc_summary %>%
    filter(QC_status == "Not_evaluated") %>%
    select(all_of(marker_key))

  excluded_X <- bind_rows(
    excluded_X,
    untested_X
  ) %>%
    distinct()

  qc_status <- paste0(
    qc_status,
    "_untested_excluded"
  )
}

# Remove excluded X markers

X <- anti_join(
  X,
  excluded_X,
  by = marker_key
)

N_X_excluded <- N_X_original - nrow(X)

# ============================================================
# 6. MATERNAL INFORMATIVE SITE SELECTION
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
      ratio_Other <= 0.10,
      ratio_Ref_mother >= 0.9 |
      ratio_Ref_mother <= 0.1
    ) %>%
    mutate(

      fetal_reads = ifelse(
        ratio_Ref_mother >= 0.9,
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
# 7. Y CHROMOSOME READS
# ============================================================

info_Y <- Y %>%
  mutate(
    Y_reads = N_REF + N_ALT,
    Covered = Y_reads > 0
  )

# ============================================================
# 8. CALCULATE MFF
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

MFF_Y <- mean(info_Y$Y_reads)

# ============================================================
# 9. CALCULATE A, B AND C RATIOS
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
# 10. FETAL SEX CLASSIFICATION
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
# 11. SUMMARY
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

  N_X_original = N_X_original,

  N_X_excluded = N_X_excluded,

  N_X_retained = nrow(X),

  N_X_sites = nrow(info_X),

  N_Y_sites = nrow(info_Y),

  N_Y_covered_sites = sum(info_Y$Covered),

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

  X_specificity_QC = qc_status,

  Min_fetal_reads = if (
    is.null(opt$min_fetal_reads)
  ) NA_integer_ else opt$min_fetal_reads,

  Min_fetal_percent = if (
    is.null(opt$min_fetal_percent)
  ) NA_real_ else opt$min_fetal_percent,

  Exclude_untested_X = opt$exclude_untested_X,

  Min_male_depth = opt$min_male_depth,

  Max_minor_fraction = opt$max_minor_fraction,

  Min_males_failed = opt$min_males_failed,

  A_threshold = opt$A_threshold,

  B_female_max = opt$B_female_max,

  B_male_min = opt$B_male_min,

  C_female_max = opt$C_female_max,

  C_male_min = opt$C_male_min
)

# ============================================================
# 12. SAVE OUTPUT FILES
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

write_output(
  excluded_X,
  "X_excluded_markers.txt"
)

if (!is.null(qc_summary)) {

  write_output(
    qc_summary,
    "X_specificity_QC.txt"
  )
}

# ============================================================
# 13. TERMINAL OUTPUT
# ============================================================

cat("\n====================================\n")
cat("NEWPAT - FETAL SEX DETERMINATION\n")
cat("====================================\n\n")

cat("Sample:", summary_result$Sample, "\n")

cat("X specificity QC:", qc_status, "\n")

cat("X markers excluded:", N_X_excluded, "\n\n")

cat("Autosomal informative sites:",
    nrow(info_A), "\n")

cat("X informative sites:",
    nrow(info_X), "\n")

cat("Y covered sites:",
    sum(info_Y$Covered), "/", nrow(info_Y), "\n\n")

cat("MFF_A:", MFF_A, "\n")
cat("MFF_X:", MFF_X, "\n")
cat("MFF_Y:", MFF_Y, "\n\n")

cat("Ratio A (X/autosomes):", A, "\n")
cat("Ratio B (Y/autosomes):", B, "\n")
cat("Ratio C (Y/X):", C, "\n\n")

cat("Fetal sex:", fetal_sex, "\n")

cat("\nResults saved in:", opt$out, "\n")

