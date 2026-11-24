#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
})

parser <- ArgumentParser(description="Postnatal CPI from SNPcall or DNPcall trio outputs")
parser$add_argument("--pileup_dir", required=TRUE)
parser$add_argument("--out_dir", required=TRUE)
parser$add_argument("--trios", required=TRUE,
                    help="mother:father:child identifiers, separated by commas")
parser$add_argument("--mother_prefix", default="full_output_")
parser$add_argument("--father_prefix", default="full_output_")
parser$add_argument("--child_prefix", default="full_output_")
parser$add_argument("--popfreq", required=TRUE)
parser$add_argument("--markers", required=FALSE)
parser$add_argument("--min_coverage", type="double", default=NULL,
                    help="Minimum total coverage per site in each trio member (REF + ALT + Other)")

parser$add_argument("--mismatch_probability_snp", default="0.001")
parser$add_argument("--mismatch_probability_dnp", default="0.0000892303778101497")
parser$add_argument("--include_hom_hom",
                    choices=c("TRUE", "FALSE"), default="TRUE")
parser$add_argument("--include_hom_het",
                    choices=c("TRUE", "FALSE"), default="TRUE")
parser$add_argument("--include_het_hom",
                    choices=c("TRUE", "FALSE"), default="TRUE")
parser$add_argument("--include_het_het",
                    choices=c("TRUE", "FALSE"), default="TRUE")

args <- parser$parse_args()

selection <- c(
  hom_hom=args$include_hom_hom == "TRUE",
  hom_het=args$include_hom_het == "TRUE",
  het_hom=args$include_het_hom == "TRUE",
  het_het=args$include_het_het == "TRUE"
)

if (!any(selection)) {
  stop("At least one genotype combination must be enabled")
}

message("Selected combinations: ",
        paste(names(selection)[selection], collapse=", "))

min_coverage <- args$min_coverage

if (!is.null(min_coverage)) {
  if (!is.finite(min_coverage) || min_coverage < 1 ||
      min_coverage != floor(min_coverage) ||
      min_coverage >= .Machine$integer.max) {
    stop("--min_coverage must be a positive integer below the R integer limit")
  }
  min_coverage <- as.integer(min_coverage)
}

message("Minimum coverage per trio member: ",
        if (is.null(min_coverage)) "not set" else min_coverage)
if (dir.exists(args$out_dir) && length(list.files(args$out_dir, all.files=TRUE,
                                                no..=TRUE))) {
  stop("Output directory is not empty; choose a new --out_dir")
}
dir.create(args$out_dir, recursive=TRUE, showWarnings=FALSE)

trio_parts <- strsplit(args$trios, ",", fixed=TRUE)[[1]]
trios <- lapply(trio_parts, function(z) trimws(strsplit(z, ":", fixed=TRUE)[[1]]))
if (!length(trios) || any(lengths(trios) != 3L) ||
    any(!nzchar(unlist(trios))) || anyDuplicated(trio_parts)) {
  stop("--trios must be mother:father:child entries separated by commas")
}
if (any(!grepl("^[A-Za-z0-9_.-]+$", unlist(trios)))) {
  stop("Identifiers may contain letters, digits, underscores, dots and hyphens")
}
if (any(vapply(trios, function(z) anyDuplicated(z) > 0, logical(1)))) {
  stop("Mother, father and child identifiers must be distinct in each trio")
}

## Read coordinates from the frequency table and use them in every join.
if (!file.exists(args$popfreq)) stop("popfreq file not found")
population <- read.table(args$popfreq, header=TRUE, sep="\t",
                         stringsAsFactors=FALSE, check.names=FALSE,
                         quote="", comment.char="")

marker_type <- function(x, label) {
  has_snp <- "pos" %in% names(x)
  has_dnp <- all(c("pos1", "pos2") %in% names(x))
  if (has_snp == has_dnp) {
    stop(label, ": use either pos for SNPs or pos1 and pos2 for DNPs")
  }
  if (has_snp) "SNP" else "DNP"
}

has_markers <- !is.null(args$markers) && nzchar(trimws(args$markers))
marker_list <- NULL
marker <- marker_type(population, "population frequencies")
if (has_markers) {
  # DNPcall lists have five columns without a header; SNPcall lists have four.
  # A header is also accepted, with lower- or upper-case allele column names.
  marker_list <- read.table(args$markers, header=FALSE, sep="", quote="",
                            comment.char="", colClasses="character",
                            stringsAsFactors=FALSE, check.names=FALSE)
  marker_columns <- switch(as.character(ncol(marker_list)),
    "4" = c("chr", "pos", "REF", "ALT"),
    "5" = c("chr", "pos1", "pos2", "REF", "ALT"),
    stop("Marker list must contain four SNP columns or five DNP columns")
  )
  first_row <- unlist(marker_list[1, ], use.names=FALSE)
  if (identical(tolower(first_row), tolower(marker_columns))) {
    marker_list <- marker_list[-1, , drop=FALSE]
  }
  names(marker_list) <- marker_columns
  if (!nrow(marker_list)) stop("Marker list is empty")
  marker <- marker_type(marker_list, "marker list")
  if (marker_type(population, "population frequencies") != marker) {
    stop("Population frequency coordinates do not match the marker list format")
  }
}
message("Panel: ", if (has_markers) args$markers else "all calling-output sites")
site_columns <- if (marker == "SNP") c("chr", "pos") else c("chr", "pos1", "pos2")
message("Marker type: ", marker)

# Select the constant for the detected marker type.
read_error_constant <- function(value, option) {
  if (is.null(value)) return(NULL)
  number <- suppressWarnings(as.numeric(value))
  if (length(number) != 1L || !is.finite(number) || number <= 0 || number > 1) {
    stop(option, " must be a finite number in (0, 1], for example 1e-04")
  }
  number
}

error_constants <- list(
  SNP = read_error_constant(args$mismatch_probability_snp, "--mismatch_probability_snp"),
  DNP = read_error_constant(args$mismatch_probability_dnp, "--mismatch_probability_dnp")
)
mismatch_probability <- error_constants[[marker]]
message("Mismatch probability used for ", marker, ": ", sprintf("%.17g", mismatch_probability))

site_key <- function(x) {
  do.call(paste, c(unname(as.list(x[site_columns])), sep="\t"))
}

check_sites <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, ": missing columns ", paste(missing, collapse=", "))
  if (anyNA(x[site_columns]) || anyDuplicated(site_key(x))) {
    stop(label, ": missing or duplicate coordinates")
  }
}

if (has_markers) {
  check_sites(marker_list, marker_columns, "marker list")
  for (column in setdiff(site_columns, "chr")) {
    position <- suppressWarnings(as.numeric(marker_list[[column]]))
    if (any(!is.finite(position) | position < 1 | position != floor(position))) {
      stop("Marker list positions must be positive integers")
    }
    marker_list[[column]] <- position
  }
  marker_list$REF <- toupper(marker_list$REF)
  marker_list$ALT <- toupper(marker_list$ALT)
  allele_pattern <- if (marker == "SNP") "^[ACGT]$" else "^[ACGT]{2}$"
  if (anyNA(marker_list[c("REF", "ALT")]) ||
      any(!grepl(allele_pattern, marker_list$REF) |
          !grepl(allele_pattern, marker_list$ALT) |
          marker_list$REF == marker_list$ALT)) {
    stop("Marker list must contain distinct REF and ALT alleles for the selected marker type")
  }
  # Check again after converting positions, so equivalent coordinates cannot repeat.
  check_sites(marker_list, marker_columns, "marker list")
  
}

check_sites(population, c(site_columns, "Ref_f.pop", "Alt_f.pop"),
            "population frequencies")
population$chr <- as.character(population$chr)
if (has_markers) {
  frequency_index <- match(site_key(marker_list), site_key(population))
  if (anyNA(frequency_index)) {
    missing_sites <- site_key(marker_list)[is.na(frequency_index)]
    stop("Missing population frequencies for marker-list sites: ",
         paste(head(missing_sites, 5), collapse=", "))
  }
  population <- population[frequency_index, , drop=FALSE]
  if (all(c("REF", "ALT") %in% names(population))) {
    if (anyNA(population[c("REF", "ALT")]) ||
        any(toupper(population$REF) != marker_list$REF |
            toupper(population$ALT) != marker_list$ALT)) {
      stop("Population frequency alleles differ from the marker list")
    }
  }
}
population$Ref_f.pop <- as.numeric(gsub(",", ".", population$Ref_f.pop))
population$Alt_f.pop <- as.numeric(gsub(",", ".", population$Alt_f.pop))
if (any(!is.finite(population$Ref_f.pop) | !is.finite(population$Alt_f.pop) |
        population$Ref_f.pop <= 0 | population$Ref_f.pop > 1 |
        population$Alt_f.pop <= 0 | population$Alt_f.pop > 1)) {
  stop("Population frequencies must be finite and in (0, 1]")
}

# Keep allele columns outside the join to avoid REF and ALT suffix conflicts.
population_for_join <- population[, c(site_columns, "Ref_f.pop", "Alt_f.pop"), drop=FALSE]

read_calls <- function(path) {
  x <- read.table(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  if (marker_type(x, path) != marker) stop(path, ": format differs from population frequencies")
  check_sites(x, c(site_columns, "REF", "ALT", "N_REF", "N_ALT", "N_Other", "total_count"), path)
  x$chr <- as.character(x$chr)
  if (has_markers) {
    # Restrict all downstream calculations to the supplied marker list.
    x <- x[site_key(x) %in% site_key(marker_list), , drop=FALSE]
    j <- match(site_key(x), site_key(marker_list))
    if (anyNA(x[c("REF", "ALT")]) ||
        any(toupper(x$REF) != marker_list$REF[j] |
            toupper(x$ALT) != marker_list$ALT[j])) {
      stop(path, ": alleles differ from the marker list")
    }
  }
  expected_length <- if (marker == "SNP") 1L else 2L
  if (anyNA(x[c("REF", "ALT")]) ||
      any(nchar(x$REF) != expected_length | nchar(x$ALT) != expected_length |
          toupper(x$REF) == toupper(x$ALT))) {
    stop(path, ": invalid REF or ALT alleles")
  }
  # Older DNP frequency tables have no REF or ALT columns.
  if (all(c("REF", "ALT") %in% names(population))) {
    j <- match(site_key(x), site_key(population))
    present <- which(!is.na(j))
    if (anyNA(population[j[present], c("REF", "ALT"), drop=FALSE])) {
      stop("Population frequencies contain missing alleles")
    }
    if (any(toupper(x$REF[present]) != toupper(population$REF[j[present]]) |
            toupper(x$ALT[present]) != toupper(population$ALT[j[present]]))) {
      stop(path, ": alleles differ from population frequencies")
    }
  }
  x
}

# Genotype intervals follow the cfDNA script. Intermediate ratios are unassigned.
classify_genotype <- function(ratio) {
  gt <- rep(NA_character_, length(ratio))
  gt[is.finite(ratio) & ratio >= 0.9] <- "0/0"
  gt[is.finite(ratio) & ratio <= 0.1] <- "1/1"
  gt[is.finite(ratio) & ratio >= 0.4 & ratio <= 0.6] <- "0/1"
  gt
}

# Apply the same Other/total filter separately to every sample.
passes_quality <- function(ref, alt, other, total, coverage=NULL) {
  ratio <- other / total

  valid <- total > 0 & (ref + alt) > 0 &
    is.finite(ratio) & ratio <= 0.10

  if (!is.null(coverage)) {
    valid <- valid & total >= coverage
  }

  valid
}
stopifnot(identical(passes_quality(c(90,89,0), c(0,0,0), c(10,11,0),
                                  c(100,100,0)), c(TRUE,FALSE,FALSE)))

classify_combination <- function(m, child) {
  case_when(
    is.na(m) | is.na(child) ~ NA_character_,
    m == "0/1" & child == "0/1" ~ "het_het",
    m == "0/1" & child %in% c("0/0", "1/1") ~ "het_hom",
    m %in% c("0/0", "1/1") & child == "0/1" ~ "hom_het",
    m == child & m %in% c("0/0", "1/1") ~ "hom_hom",
    TRUE ~ NA_character_
  )
}

required_allele <- function(m, child) {
  case_when(
    m == "0/0" & child == "0/1" ~ "ALT",
    m == "1/1" & child == "0/1" ~ "REF",
    m %in% c("0/0", "0/1") & child == "0/0" ~ "REF",
    m %in% c("1/1", "0/1") & child == "1/1" ~ "ALT",
    TRUE ~ NA_character_
  )
}

prepare_sample <- function(id, prefix, role) {
  path <- file.path(args$pileup_dir, paste0(prefix, id, ".txt"))
  x <- read_calls(path)
  pattern <- if (marker == "SNP") "^[ACGT]$" else "^[ACGT]{2}$"
  if (any(!grepl(pattern, toupper(x$REF)) | !grepl(pattern, toupper(x$ALT)))) {
    stop(path, ": alleles must contain only A, C, G and T")
  }
  count_columns <- c("N_REF", "N_ALT", "N_Other", "total_count")
  if (any(!vapply(x[count_columns], is.numeric, logical(1)))) {
    stop(path, ": read counts must be numeric")
  }
  # DNPcall can report NA total_count at sites with no reads.
  # Correct this representation in memory and retain an internal flag.
  x$Missing_total_at_zero_coverage <- is.na(x$total_count) &
    !is.na(x$N_REF) & !is.na(x$N_ALT) & !is.na(x$N_Other) &
    x$N_REF == 0 & x$N_ALT == 0 & x$N_Other == 0
  if (any(x$Missing_total_at_zero_coverage)) {
    message(path, ": ", sum(x$Missing_total_at_zero_coverage),
            " zero-coverage site(s) with missing total_count; treated as zero coverage")
    x$total_count[x$Missing_total_at_zero_coverage] <- 0
  }
  counts <- as.matrix(x[count_columns])
  if (anyNA(counts) || any(!is.finite(counts) | counts < 0 | counts != floor(counts))) {
    stop(path, ": counts must be nonnegative finite integers")
  }
  if (any(x$total_count != x$N_REF + x$N_ALT + x$N_Other)) {
    stop(path, ": total_count differs from REF + ALT + Other")
  }
  x$ratio <- x$N_REF / (x$N_REF + x$N_ALT)
  x$GT <- classify_genotype(x$ratio)
  x$Other_ratio <- x$N_Other / x$total_count
  x$Quality_pass <- passes_quality(
    x$N_REF, x$N_ALT, x$N_Other, x$total_count,
    coverage=min_coverage
  )
  x$GT[!x$Quality_pass] <- NA_character_
  x$Present <- TRUE
  x <- x[c(site_columns, "REF", "ALT", count_columns, "ratio", "GT",
           "Other_ratio", "Quality_pass", "Present", "Missing_total_at_zero_coverage")]
  names(x)[!names(x) %in% site_columns] <- paste0(
    names(x)[!names(x) %in% site_columns], "_", role)
  x
}

# PI for sites with an obligatory paternal allele.
# The same formula applies to compatible homozygous
# and heterozygous mothers.

site_pi <- function(required, father_gt, ref_frequency, alt_frequency, epsilon) {
  frequency <- ifelse(required == "REF", ref_frequency, alt_frequency)
  transmission <- ifelse(father_gt == "0/1", 0.5,
    ifelse((required == "REF" & father_gt == "0/0") |
           (required == "ALT" & father_gt == "1/1"), 1, 0))
  mismatch <- transmission == 0
  weight <- ifelse(mismatch, epsilon, transmission)
  list(PI=weight / frequency, Mismatch=mismatch,
       Frequency=frequency, Transmission=transmission)
}

# Basic genotype and formula checks.
stopifnot(identical(classify_genotype(c(1, 0.9, 0.6, 0.5, 0.4, 0.1, 0, 0.3, NaN)),
                    c("0/0", "0/0", "0/1", "0/1", "0/1", "1/1", "1/1", NA_character_, NA_character_)))
# Include homozygous children and paternal mismatches.
test <- site_pi(c("ALT", "ALT", "REF", "REF"),
                c("1/1", "0/1", "0/0", "1/1"),
                rep(0.6, 4), rep(0.4, 4), 0.001)
stopifnot(isTRUE(all.equal(test$PI, c(2.5, 1.25, 1/0.6, 0.001/0.6))),
          identical(test$Mismatch, c(FALSE, FALSE, FALSE, TRUE)))

summary_rows <- list()

for (ids in trios) {
  trio_name <- paste0("Mother", ids[1], "-Father", ids[2], "-Child", ids[3],
                      "_", marker, "s")
  message("Processing ", trio_name)
  mother <- prepare_sample(ids[1], args$mother_prefix, "mother")
  father <- prepare_sample(ids[2], args$father_prefix, "father")
  child <- prepare_sample(ids[3], args$child_prefix, "child")
  x <- full_join(full_join(mother, father, by=site_columns), child, by=site_columns)
  if (has_markers) {
    x <- full_join(marker_list[site_columns], x, by=site_columns)
  }

  # Check all shared allele definitions, including sites without frequencies.
  for (pair in list(c("mother", "father"), c("mother", "child"), c("father", "child"))) {
    for (allele in c("REF", "ALT")) {
      a <- x[[paste0(allele, "_", pair[1])]]
      b <- x[[paste0(allele, "_", pair[2])]]
      if (any(!is.na(a) & !is.na(b) & toupper(a) != toupper(b))) {
        stop(trio_name, ": inconsistent allele definitions between samples")
      }
    }
  }
  x <- left_join(x, population_for_join, by=site_columns)
  x$REF <- coalesce(x$REF_mother, x$REF_father, x$REF_child)
  x$ALT <- coalesce(x$ALT_mother, x$ALT_father, x$ALT_child)
  if (has_markers) {
    panel_index <- match(site_key(x), site_key(marker_list))
    x$REF <- coalesce(x$REF, marker_list$REF[panel_index])
    x$ALT <- coalesce(x$ALT, marker_list$ALT[panel_index])
  }
  m <- x$GT_mother
  cgt <- x$GT_child

  x$Mother_child_incompatible <- !is.na(m) & !is.na(cgt) &
    ((m == "0/0" & cgt == "1/1") |
     (m == "1/1" & cgt == "0/0"))

  x$Genotype_combination <- classify_combination(m, cgt)

  x$Required_paternal_allele <- required_allele(m, cgt)

  x$Combination_selected <-
    x$Genotype_combination %in% names(selection)[selection]

  x$Status <- case_when(
    is.na(x$Present_mother) |
      is.na(x$Present_father) |
      is.na(x$Present_child) ~ "missing_call",

    !x$Quality_pass_mother |
      !x$Quality_pass_father |
      !x$Quality_pass_child ~ "quality_filter",

    is.na(m) ~ "mother_unassigned",
    is.na(cgt) ~ "child_unassigned",
    is.na(x$GT_father) ~ "father_unassigned",

    x$Mother_child_incompatible ~ "mother_child_incompatible",

    !x$Combination_selected ~ "combination_not_selected",

    is.na(x$Ref_f.pop) |
      is.na(x$Alt_f.pop) ~ "missing_population_frequency",

    TRUE ~ "eligible"
  )

  eligible <- x$Status == "eligible"

  obligatory <- eligible &
    !is.na(x$Required_paternal_allele)

  no_obligatory <- eligible &
    x$Genotype_combination %in% "het_het"

  stopifnot(all(!eligible | obligatory | no_obligatory))

  x$Paternity_Index <- NA_real_
  x$Mismatch <- NA
  x$Paternal_transmission <- NA_real_
  x$Paternity_status <- NA_character_

  if (any(obligatory)) {
    pi <- site_pi(
      x$Required_paternal_allele[obligatory],
      x$GT_father[obligatory],
      x$Ref_f.pop[obligatory],
      x$Alt_f.pop[obligatory],
      mismatch_probability
    )

    x$Paternity_Index[obligatory] <- pi$PI
    x$Mismatch[obligatory] <- pi$Mismatch
    x$Paternal_transmission[obligatory] <- pi$Transmission

    x$Paternity_status[obligatory] <-
      ifelse(pi$Mismatch, "mismatch", "match")
  }

  # Mother and child are both heterozygous.
  # PI equals 1 when REF and ALT frequencies sum to 1.
  if (any(no_obligatory)) {
    x$Paternity_Index[no_obligatory] <-
      1 / (x$Ref_f.pop[no_obligatory] +
           x$Alt_f.pop[no_obligatory])

    x$Paternity_status[no_obligatory] <-
      "no_obligatory_allele"
  }

  sites <- x[eligible, , drop=FALSE]

  log_cpi <- if (nrow(sites)) {
    sum(log10(sites$Paternity_Index))
  } else {
    NA_real_
  }

  cpi <- if (is.na(log_cpi)) NA_real_ else 10^log_cpi

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    Trio=trio_name,
    Marker_type=marker,

    N_sites=nrow(sites),

    N_hom_hom=sum(sites$Genotype_combination == "hom_hom"),
    N_hom_het=sum(sites$Genotype_combination == "hom_het"),
    N_het_hom=sum(sites$Genotype_combination == "het_hom"),
    N_het_het=sum(sites$Genotype_combination == "het_het"),

    Matching_sites=sum(sites$Mismatch %in% FALSE),
    Mismatching_sites=sum(sites$Mismatch %in% TRUE),

    No_obligatory_allele_sites=sum(no_obligatory),

    Maternal_incompatible_sites=
      sum(x$Status == "mother_child_incompatible"),

    Log10_CPI=log_cpi,
    CPI=cpi,

    Mismatch_probability=mismatch_probability,

    Min_coverage=if (is.null(min_coverage)) {
      NA_integer_
    } else {
      min_coverage
    },

    Include_hom_hom=unname(selection["hom_hom"]),
    Include_hom_het=unname(selection["hom_het"]),
    Include_het_hom=unname(selection["het_hom"]),
    Include_het_het=unname(selection["het_het"])
  )

  write.table(
    x,
    file.path(args$out_dir, paste0(trio_name, "_sites.txt")),
    sep="\t", quote=FALSE, row.names=FALSE, na="NA"
  )
}

summary <- bind_rows(summary_rows)
write.table(summary, file.path(args$out_dir, "Summary_Postnatal_CPI.txt"),
            sep="\t", quote=FALSE, row.names=FALSE, na="NA")

message("Done: trio CPI summary and site-level results saved.")
