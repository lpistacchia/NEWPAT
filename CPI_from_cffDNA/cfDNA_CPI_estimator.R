#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(argparse)
  library(stringr)
  library(tidyr)
})

## cli parser
parser <- ArgumentParser()

parser$add_argument("--pileup_dir",    required=TRUE)
parser$add_argument("--out_dir",       required=TRUE)
parser$add_argument("--pairs",         required=TRUE)
parser$add_argument("--mother_prefix", required=TRUE)
parser$add_argument("--father_prefix", required=TRUE)
parser$add_argument("--popfreq",       required=TRUE)

## optional cutoff and err_const
parser$add_argument("--cutoff_map",    required=FALSE,
                    help='optional mapping "0.3=/abs/fileA.xlsx,0.5=/abs/fileB.xlsx". if missing use ALL sites')
parser$add_argument("--err_const",     default="0.0000892303778101497")

args <- parser$parse_args()
err_const <- as.numeric(args$err_const)

## make output dirs
if (!dir.exists(args$out_dir)) dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(args$out_dir, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

## parse pairs
pairs_vec <- strsplit(args$pairs, ",")[[1]] |> trimws()
coppie <- do.call(rbind, lapply(pairs_vec, function(x){
  sp <- strsplit(x, ":", fixed=TRUE)[[1]] |> trimws()
  data.frame(mamma=sp[1], papa=sp[2])
}))
rownames(coppie) <- NULL

## parse cutoff map OR set to ALL mode
if (!is.null(args$cutoff_map) && nzchar(args$cutoff_map)) {

  items <- strsplit(args$cutoff_map, ",")[[1]] |> trimws()

  cutoff_list <- lapply(items, function(x){
    if (grepl("=", x)) {
      ## case: label=path
      kv <- strsplit(x, "=", fixed=TRUE)[[1]] |> trimws()
      if (length(kv) != 2) stop(paste0("invalid cutoff_map entry: ", x))
      list(Cutoff = kv[1], Path = kv[2])
    } else {
      ## case: only filename provided → use filename as label too
      list(Cutoff = x, Path = x)
    }
  })

  cutoff_df <- do.call(rbind, lapply(cutoff_list, as.data.frame))
  has_cutoffs <- TRUE

} else {
  cutoff_df <- data.frame(Cutoff="ALL", Path=NA_character_)
  has_cutoffs <- FALSE
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

marker <- marker_type(population, "population frequencies")
site_columns <- if (marker == "SNP") c("chr", "pos") else c("chr", "pos1", "pos2")
message("Marker type: ", marker)

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

check_sites(population, c(site_columns, "Ref_f.pop", "Alt_f.pop"),
            "population frequencies")
population$chr <- as.character(population$chr)
population$Ref_f.pop <- as.numeric(gsub(",", ".", population$Ref_f.pop))
population$Alt_f.pop <- as.numeric(gsub(",", ".", population$Alt_f.pop))
if (any(!is.finite(population$Ref_f.pop) | !is.finite(population$Alt_f.pop) |
        population$Ref_f.pop <= 0 | population$Ref_f.pop > 1 |
        population$Alt_f.pop <= 0 | population$Alt_f.pop > 1)) {
  stop("Population frequencies must be finite and in (0, 1]")
}

# Keep allele columns outside the join to avoid REF and ALT suffix conflicts.
PopFreqDNPs <- population[, c(site_columns, "Ref_f.pop", "Alt_f.pop"), drop=FALSE]

read_calls <- function(path) {
  x <- read.table(path, header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  if (marker_type(x, path) != marker) stop(path, ": format differs from population frequencies")
  check_sites(x, c(site_columns, "REF", "ALT", "N_REF", "N_ALT", "N_Other", "total_count"), path)
  x$chr <- as.character(x$chr)
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

check_pair <- function(mother, father) {
  j <- match(site_key(mother), site_key(father))
  present <- which(!is.na(j))
  if (any(toupper(mother$REF[present]) != toupper(father$REF[j[present]]) |
          toupper(mother$ALT[present]) != toupper(father$ALT[j[present]]))) {
    stop("Mother and father have different REF or ALT alleles at shared coordinates")
  }
}

read_cutoff <- function(path) {
  if (is.na(path)) return(NULL)
  x <- read.table(path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if (marker_type(x, path) != marker) stop(path, ": cutoff format differs from the calls")
  check_sites(x, site_columns, path)
  x$chr <- as.character(x$chr)
  x
}

## helper: compute fetal_reads on mother dataframe
## note: this mirrors the logic used in cpi iteration
compute_madre_info <- function(madre_df) {
  madre_df %>%
    dplyr::filter(
      (ratio_Ref_mother >= 0.9 & Alt_mother >= 2) |
        (ratio_Ref_mother <= 0.1 & Ref_mother >= 2)
    ) %>%
    dplyr::mutate(
      fetal_reads = ifelse(ratio_Ref_mother >= 0.9, Alt_mother, Ref_mother)
    )
}

## core: compute cpi for one list (or all sites if dnp_list is NULL)
run_CPI_iteration <- function(madre_df, padre_df, dnp_list, cutoff_label, err_const) {
  if (!is.null(dnp_list)) {
    madre <- madre_df %>% dplyr::semi_join(dnp_list, by=site_columns)
    padre <- padre_df %>% dplyr::semi_join(dnp_list, by=site_columns)
  } else {
    madre <- madre_df
    padre <- padre_df
  }

  madre_info <- compute_madre_info(madre)
  max_fetal <- suppressWarnings(max(madre_info$fetal_reads, na.rm=TRUE))
  if (!is.finite(max_fetal) || max_fetal < 2)
    return(list(results=data.frame(), madre_info=madre_info))

  results <- data.frame()
  for (i in 2:max_fetal) {
    madre_filt <- madre_info %>% dplyr::filter(fetal_reads >= i)

    df_merged <- dplyr::inner_join(madre_filt, padre, by=site_columns) %>%
      dplyr::mutate(threshold = dplyr::case_when(
        ratio_Ref_mother >= 0.9 & Alt_father >= 0 ~ "informative",
        ratio_Ref_mother <= 0.1 & Ref_father >= 0 ~ "informative",
        TRUE ~ "non informative"
      )) %>%
      dplyr::inner_join(PopFreqDNPs, by=site_columns) %>%
      dplyr::mutate(dplyr::across(c(Ref_f.pop,Alt_f.pop), as.numeric))

    PI <- df_merged %>% dplyr::filter(threshold=="informative") %>%
      dplyr::mutate(Paternity_Index = dplyr::case_when(
        ratio_Ref_mother >= 0.9 & ratio_Ref_father <= 0.1 ~ 1/Alt_f.pop,
        ratio_Ref_mother >= 0.9 & between(ratio_Ref_father,0.4,0.6) ~ 0.5/Alt_f.pop,
        ratio_Ref_mother >= 0.9 & ratio_Ref_father >= 0.9 ~ err_const/Alt_f.pop,
        ratio_Ref_mother <= 0.1 & ratio_Ref_father <= 0.1 ~ err_const/Ref_f.pop,
        ratio_Ref_mother <= 0.1 & between(ratio_Ref_father,0.4,0.6) ~ 0.5/Ref_f.pop,
        ratio_Ref_mother <= 0.1 & ratio_Ref_father >= 0.9 ~ 1/Ref_f.pop,
        TRUE ~ NA_real_
      )) %>% tidyr::drop_na()

    CPI_val <- ifelse(nrow(PI)>0, prod(PI$Paternity_Index, na.rm=TRUE), NA_real_)
    logCPI  <- ifelse(!is.na(CPI_val) & CPI_val>0, log10(CPI_val), NA_real_)

    mismatch_summary <- df_merged %>% dplyr::filter(threshold=="informative") %>%
      summarise(
        Mismatch = sum((ratio_Ref_mother>=0.9 & ratio_Ref_father>=0.9) |
                         (ratio_Ref_mother<=0.1 & ratio_Ref_father<=0.1)),
        NotMismatch = n() - Mismatch
      )

    results <- rbind(results, data.frame(
      Iteration=i,
      InfSites_mother=nrow(madre_filt),
      InfSites_father=nrow(PI),
      CPI=CPI_val,
      log_CPI=logCPI,
      Mismatch=mismatch_summary$Mismatch,
      NotMismatch=mismatch_summary$NotMismatch,
      Cutoff=cutoff_label
    ))
  }
  return(list(results=results, madre_info=madre_info))
}

## get logCPI at iteration closest to target (median fetal reads)
closest_logCPI <- function(res_df, target_i) {
  if (nrow(res_df)==0 || is.na(target_i) || !is.finite(target_i)) return(NA_real_)
  idx <- which.min(abs(res_df$Iteration - target_i))
  res_df$log_CPI[idx]
}

plot_CPI <- function(df, cutoff_label, pair_name,
                     cfDNA_percent = NA_real_,
                     cpi_medff = NA_real_) {

  if (nrow(df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste(pair_name, "-", cutoff_label, "(no data)"))
    )
  }

  max_iter <- max(df$Iteration, na.rm = TRUE)

  ticks <- if (max_iter <= 20) {
    2:max_iter
  } else if (max_iter <= 100) {
    unique(c(2, seq(10, max_iter, by = 10), max_iter))
  } else if (max_iter <= 500) {
    unique(c(2, seq(25, max_iter, by = 25), max_iter))
  } else {
    unique(c(2, seq(50, max_iter, by = 50), max_iter))
  }

  cfDNA_lbl <- if (is.finite(cfDNA_percent)) {
    sprintf("%.2f", cfDNA_percent)
  } else {
    "NA"
  }

  cpi_lbl <- if (is.finite(cpi_medff)) {
    sprintf("%.3f", cpi_medff)
  } else {
    "NA"
  }

  annot_text <- paste0(
    "atop(",
    "'cfDNA: ", cfDNA_lbl, "%',",
    "log[10](CPI) * ': ' * '", cpi_lbl, "'",
    ")"
  )

  ggplot(df, aes(x = Iteration)) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = -4, fill = "red", alpha = 0.05) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -4, ymax = 4, fill = "grey", alpha = 0.05) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = 4, ymax = Inf, fill = "green", alpha = 0.05) +

    geom_line(aes(y = log_CPI), colour = "black") +
    geom_point(aes(y = log_CPI), colour = "black", size = 1.3) +
    geom_line(aes(y = NotMismatch), colour = "green") +
    geom_line(aes(y = Mismatch), colour = "orange",
              linetype = "dashed") +

    annotate(
      "label",
      x = Inf, y = Inf,
      label = annot_text,
      parse = TRUE,
      hjust = 1.05,
      vjust = 1.10,
      size = 3,
      fill = "white",
      colour = "black",
      label.size = 0.2
    ) +

    scale_y_continuous(
      name = expression(log[10](CPI)),
      sec.axis = sec_axis(~ ., name = "Matching / mismatching sites"),
      expand = expansion(mult = c(0.02, 0.14))
    ) +

    scale_x_continuous(
      name = "minimum number of fetal reads considered",
      breaks = ticks,
      limits = c(2, max(3, max_iter)),
      expand = expansion(mult = c(0.02, 0.02))
    ) +

    labs(title = paste(pair_name, "-", cutoff_label)) +

    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(size = 9, hjust = 0),
      axis.title = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.6
      ),
      plot.margin = margin(8, 12, 8, 8)
    )
}

## main containers
summary_table <- data.frame()
all_combined_plots <- list()

for (i in 1:nrow(coppie)) {
  mamma_id <- coppie$mamma[i]
  papa_id  <- coppie$papa[i]
  pair_name <- paste0("Mother", mamma_id, "-Father", papa_id, "_", marker, "s")

  mother_path <- file.path(args$pileup_dir, paste0(args$mother_prefix, mamma_id, ".txt"))
  father_path <- file.path(args$pileup_dir, paste0(args$father_prefix, papa_id, ".txt"))

  madre <- read_calls(mother_path)
  padre <- read_calls(father_path)
  check_pair(madre, padre)

  madre$chr <- as.character(madre$chr)
  padre$chr <- as.character(padre$chr)

  madre <- madre %>%
    dplyr::transmute(
      dplyr::across(dplyr::all_of(site_columns)), REF, ALT,
      Ref_mother = N_REF,
      Alt_mother = N_ALT,
      Error_mother = N_Other,
      ratio_Ref_mother = N_REF / (N_REF + N_ALT),
      total_count = total_count
    )

  padre <- padre %>%
    dplyr::transmute(
      dplyr::across(dplyr::all_of(site_columns)), REF, ALT,
      Ref_father = N_REF,
      Alt_father = N_ALT,
      Error_father = N_Other,
      ratio_Ref_father = N_REF / (N_REF + N_ALT),
      total_count = total_count
    )

  madre <- madre %>% dplyr::filter(Error_mother/total_count <= 0.10)
  padre <- padre %>%
    dplyr::mutate(
      Ref_father = ifelse(ratio_Ref_father < 0.1, 0, Ref_father),
      Alt_father = ifelse(ratio_Ref_father > 0.9, 0, Alt_father)
    ) %>% dplyr::filter(Error_father/total_count <= 0.10)

  per_cutoff_results <- list()
  per_cutoff_plots   <- list()
  per_cutoff_CPI_medreads <- c()

  for (k in seq_len(nrow(cutoff_df))) {
    cutoff_label <- cutoff_df$Cutoff[k]
    dnp_list <- read_cutoff(cutoff_df$Path[k])

    out <- run_CPI_iteration(madre, padre, dnp_list, cutoff_label, err_const)
    res <- out$results
    per_cutoff_results[[cutoff_label]] <- res

    target_i <- out$madre_info %>%
      dplyr::pull(fetal_reads) %>%
      stats::median(na.rm = TRUE) %>%
      round()

    cpi_at_med <- closest_logCPI(res, target_i)

    per_cutoff_CPI_medreads <- c(
      per_cutoff_CPI_medreads,
      cpi_at_med
    )

    median_fetal_fraction <- out$madre_info %>%
      dplyr::mutate(
        fetal_fraction_site = 2 * fetal_reads /
          (Ref_mother + Alt_mother)
      ) %>%
      dplyr::pull(fetal_fraction_site) %>%
      stats::median(na.rm = TRUE)

    cfDNA_percent <- median_fetal_fraction * 100

    per_cutoff_plots[[cutoff_label]] <- plot_CPI(
      res,
      cutoff_label,
      pair_name,
      cfDNA_percent = cfDNA_percent,
      cpi_medff = cpi_at_med
    )
  }

  ## save combined plots in plots_dir
  if (length(per_cutoff_plots)>0) {
    combined <- Reduce(`+`, per_cutoff_plots[names(per_cutoff_plots)]) +
      patchwork::plot_layout(ncol=length(per_cutoff_plots))
    ggsave(file.path(plots_dir, paste0(pair_name,"_CPIplot.pdf")),
           plot=combined, width=max(12,4*length(per_cutoff_plots)), height=4)
    all_combined_plots[[pair_name]] <- combined
  }

  ## build first summary row
  all_res <- dplyr::bind_rows(per_cutoff_results)
  all_res$log_CPI <- as.numeric(all_res$log_CPI)

  summary_table <- rbind(summary_table, data.frame(
    Pair = pair_name,
    Mother_sites_max = suppressWarnings(max(all_res$InfSites_mother, na.rm=TRUE)),
    Father_sites_max = suppressWarnings(max(all_res$InfSites_father, na.rm=TRUE)),
    Max_logCPI = suppressWarnings(round(max(all_res$log_CPI, na.rm=TRUE),3)),
    Min_logCPI = suppressWarnings(round(min(all_res$log_CPI, na.rm=TRUE),3)),
    Median_logCPI = suppressWarnings(round(stats::median(all_res$log_CPI, na.rm=TRUE),3)),
    CPI_medianFF = suppressWarnings(round(stats::median(per_cutoff_CPI_medreads, na.rm=TRUE),3))
  ))

}

## save single big pdf with all pairs
if (length(all_combined_plots)>0) {
  final_plot <- wrap_plots(all_combined_plots, ncol=1)
  ggsave(file.path(plots_dir, "All_CPIplots.pdf"),
         plot=final_plot, width=12, height=4*length(all_combined_plots), limitsize=FALSE)
}


## fetal fraction + cpi summaries (cutoff or ALL)
fetal_fraction_summary <- data.frame()

for (i in 1:nrow(coppie)) {
  mamma_id <- coppie$mamma[i]
  papa_id  <- coppie$papa[i]
  pair_name <- paste0("Mother", mamma_id, "-Father", papa_id, "_", marker, "s")

  mother_path <- file.path(args$pileup_dir, paste0(args$mother_prefix, mamma_id, ".txt"))
  father_path <- file.path(args$pileup_dir, paste0(args$father_prefix, papa_id, ".txt"))

  madre <- read_calls(mother_path)
  padre <- read_calls(father_path)
  check_pair(madre, padre)

madre$chr <- as.character(madre$chr)
padre$chr <- as.character(padre$chr)

madre <- madre %>%
  dplyr::transmute(
    dplyr::across(dplyr::all_of(site_columns)), REF, ALT,
    Ref_mother = N_REF,
    Alt_mother = N_ALT,
    Error_mother = N_Other,
    ratio_Ref_mother = N_REF / (N_REF + N_ALT),
    total_count = total_count
  )

padre <- padre %>%
  dplyr::transmute(
    dplyr::across(dplyr::all_of(site_columns)), REF, ALT,
    Ref_father = N_REF,
    Alt_father = N_ALT,
    Error_father = N_Other,
    ratio_Ref_father = N_REF / (N_REF + N_ALT),
    total_count = total_count
  )

  madre <- madre %>% filter(Error_mother/total_count <= 0.10)
  padre <- padre %>%
    mutate(
      Ref_father = ifelse(ratio_Ref_father < 0.1,0,Ref_father),
      Alt_father = ifelse(ratio_Ref_father > 0.9,0,Alt_father)
    ) %>% filter(Error_father/total_count <= 0.10)

  for (k in seq_len(nrow(cutoff_df))) {
    cutoff_label <- cutoff_df$Cutoff[k]
    dnp_list <- read_cutoff(cutoff_df$Path[k])

    madre_cutoff <- if (!is.null(dnp_list)) madre %>% semi_join(dnp_list, by=site_columns) else madre
    padre_cutoff <- if (!is.null(dnp_list)) padre %>% semi_join(dnp_list, by=site_columns) else padre

    madre_cutoff <- madre_cutoff %>%
      mutate(
        fetal_reads = case_when(
          ratio_Ref_mother >= 0.9 ~ Alt_mother,
          ratio_Ref_mother <= 0.1 ~ Ref_mother,
          TRUE ~ NA_real_
        ),
        fetal_fraction_site = 2*fetal_reads/(Ref_mother+Alt_mother)
      ) %>% filter(!is.na(fetal_fraction_site), fetal_fraction_site > 0, fetal_reads >= 2)

    n_sites_mother <- nrow(madre_cutoff)

    ## recalc father informative count after merge
    if (!is.null(dnp_list)) {
      padre_cut_tmp <- padre %>% semi_join(dnp_list, by=site_columns)
    } else {
      padre_cut_tmp <- padre
    }

    df_merge_tmp <- inner_join(madre_cutoff, padre_cut_tmp, by=site_columns) %>%
      mutate(threshold = case_when(
        ratio_Ref_mother >= 0.9 & Alt_father >= 0 ~ "informative",
        ratio_Ref_mother <= 0.1 & Ref_father >= 0 ~ "informative",
        TRUE ~ "non informative"
      ))

    n_sites_father <- sum(df_merge_tmp$threshold=="informative")


    median_fetal_fraction <- suppressWarnings(median(madre_cutoff$fetal_fraction_site, na.rm=TRUE))

    out <- run_CPI_iteration(madre, padre, dnp_list, cutoff_label, err_const)
    res <- out$results
    all_cpi_vals <- res$log_CPI[!is.na(res$log_CPI)]

    max_cpi <- ifelse(length(all_cpi_vals)==0, NA, max(all_cpi_vals))
    min_cpi <- ifelse(length(all_cpi_vals)==0, NA, min(all_cpi_vals))
    med_cpi <- ifelse(length(all_cpi_vals)==0, NA, median(all_cpi_vals))
    mean_cpi <- ifelse(length(all_cpi_vals)==0, NA, mean(all_cpi_vals))

    target_i <- out$madre_info %>% pull(fetal_reads) %>% median(na.rm=TRUE) %>% round()
    cpi_at_med <- closest_logCPI(res, target_i)

    if (has_cutoffs) {
      fetal_fraction_summary <- rbind(fetal_fraction_summary, data.frame(
        Pair = pair_name,
        Cutoff = cutoff_label,
        N_sites_mother = n_sites_mother,
        N_sites_father = n_sites_father,
        Median_fetal_fraction_raw = round(median_fetal_fraction,4),
        Median_fetal_fraction_percent = round(median_fetal_fraction*100,2),
        Max_logCPI = round(max_cpi,3),
        Min_logCPI = round(min_cpi,3),
        Median_logCPI = round(med_cpi,3),
        Mean_logCPI = round(mean_cpi,3),
        CPI_medianFF = round(cpi_at_med,3)
      ))
    } else {
      fetal_fraction_summary <- rbind(fetal_fraction_summary, data.frame(
        Pair = pair_name,
        N_sites_mother = n_sites_mother,
        N_sites_father = n_sites_father,
        Median_fetal_fraction_raw = round(median_fetal_fraction,4),
        Median_fetal_fraction_percent = round(median_fetal_fraction*100,2),
        Max_logCPI = round(max_cpi,3),
        Min_logCPI = round(min_cpi,3),
        Median_logCPI = round(med_cpi,3),
        Mean_logCPI = round(mean_cpi,3),
        CPI_medianFF = round(cpi_at_med,3)
      ))
    }
  }
}

write.csv(
  fetal_fraction_summary,
  file.path(args$out_dir, "Summary_FetalFraction_CPI.csv"),
  row.names = FALSE
)

message("done: all summaries and plots saved.")
