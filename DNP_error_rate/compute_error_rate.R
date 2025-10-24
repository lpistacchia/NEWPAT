#!/usr/bin/env Rscript

# load packages quietly
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(stringr)
  library(readxl)
  library(dplyr)
})


# define command line options
option_list <- list(
  make_option(c("--pileup"), type="character", default=NULL,
              help="path to a single pileup file", metavar="character"),
  make_option(c("--pileup_dir"), type="character", default=NULL,
              help="path to a directory with multiple pileup files", metavar="character"),
  make_option(c("--dnp_list"), type="character", default="DNP_list_FINAL.txt",
              help="path to DNP list file", metavar="character"),
  make_option(c("--output"), type="character", default="./",
              help="output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# determine input files
if (!is.null(opt$pileup)) {
  files <- opt$pileup
} else if (!is.null(opt$pileup_dir)) {
  files <- list.files(path = opt$pileup_dir, pattern = "\\.txt$", full.names = TRUE)
} else {
  stop("please provide either --pileup or --pileup_dir")
}

# set output directory
output_dir <- opt$output
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# load DNP list
r_DNPlist <- read.table(opt$dnp_list, header=TRUE)
r_DNPlist$pos1 <- as.numeric(r_DNPlist$pos1)
r_DNPlist$pos2 <- as.numeric(r_DNPlist$pos2)


# function to process each pileup file
process_pileup_file <- function(file) {
  base_name <- tools::file_path_sans_ext(basename(file))

  # read pileup file
  r_bamlist <- read.csv(file, quote = "", row.names = NULL, stringsAsFactors = FALSE, header = FALSE, sep = "\t")

  # create bedfile
  bedfile <- as.data.frame(cbind(r_DNPlist$chr, r_DNPlist$pos1-2, r_DNPlist$pos2+1))

  # create dataframe to have unique ids for dnps and reads
  ForUniqReads1 <- as.data.frame(cbind(paste0(r_DNPlist$chr,":",r_DNPlist$pos1-1),paste0(r_DNPlist$chr,":",r_DNPlist$pos1),
                                       paste0(r_DNPlist$chr,":",r_DNPlist$pos2), paste0(r_DNPlist$chr,":",r_DNPlist$pos2+1),
                                       seq.int(nrow(r_DNPlist))))
  df_int <- gather(data=ForUniqReads1,key="col", value="pos", 1:4)
  ForUniqReads <- df_int[,c(3,1)]

  colnames(ForUniqReads) <- c("df.chr","uniq")

  # create ref/alt dnp dataframe
  RefAltDNP <- as.data.frame(cbind(paste0(r_DNPlist$chr,":",r_DNPlist$pos1,"|",r_DNPlist$chr,":",r_DNPlist$pos2),
                                   gsub('^(.{1})(.*)$', '\\1|\\2',r_DNPlist$REF),
                                   gsub('^(.{1})(.*)$', '\\1|\\2',r_DNPlist$ALT)))
  colnames(RefAltDNP)<-c("POS","ref","alt")


  # replace "." or "," in column V5 with the corresponding reference base from column V3
  #r_bamlist$V5 <- mapply(function(v3, v5) gsub("[.,]", v3, v5), r_bamlist$V3, r_bamlist$V5)

  # create a "chr:position" column
  r_bamlist$chr <- paste(r_bamlist$V1, r_bamlist$V2, sep = ":")

  # create df2 dataframe with the columns of interest
  df2 <- data.frame(df.chr = r_bamlist$chr, df.V5 = r_bamlist$V5, df.V7 = r_bamlist$V7)

  df2$df.V7 <- paste0(",",df2$df.V7,",")

  # replace indels with special characters
  df2$df.V5 = gsub('\\*', '@', df2$df.V5)
  df2$df.V5 = gsub('[.,actgACTG]([+-]\\d+)+', '@', df2$df.V5)

  df2$df.V5 = gsub('', ' ', df2$df.V5)
  df2$df.V7 = gsub(',', ' ', df2$df.V7)

  # split values in these columns, pairing each base with the corresponding read
  df3 <- df2 %>%
    separate_rows(df.V5, df.V7, sep = " ")

  # combine read information with base information
  df5 <- df3 %>%
    group_by(df.V7) %>%
    summarise(df.V5 = paste(df.V5, collapse = "|"))

  # combine read information with position information
  df6 <- df3 %>%
    group_by(df.V7) %>%
    summarise(df.chr = paste(df.chr, collapse = "|"))

  # merge the two datasets created in the previous steps
  df_merged <- merge(df5,df6)


  # remove from the dataset everything that is not called on both positions of a single read
  df_merged2 <- df_merged[nchar(df_merged$df.V5)>6,]
  df_merged3 <- df_merged2[nchar(df_merged2$df.V7)>0,]

  # replace characters from different strands to make them all consistent
  df_merged3$df.V5 = gsub('g', 'G', df_merged3$df.V5)
  df_merged3$df.V5 = gsub('a', 'A', df_merged3$df.V5)
  df_merged3$df.V5 = gsub('c', 'C', df_merged3$df.V5)
  df_merged3$df.V5 = gsub('t', 'T', df_merged3$df.V5)
  df_merged3$df.V5 = gsub(',', '.', df_merged3$df.V5)

  conteggio <- sum(grepl("@", df_merged3$df.V5))
  print(conteggio)

  # remove all rows where "@" is present
  df_merged4 <- df_merged3 %>%
    filter_all(all_vars(!grepl("@", .)))

  df_merged5 <- df_merged4 %>%
    mutate(df.V5.1 = substr(df.V5, 1, 2),
           dfV5.2 = substr(df.V5, 3, 5),
           dfV5.3 = substr(df.V5, 6, 7))
  df_merged5 <- subset(df_merged5, select = -df.V5)


  # count all genotypes found and remove the column with info of the 4 bases
  df_count2 <- df_merged5 %>%
    group_by(df.chr, dfV5.2) %>%
    summarise(n = n())

  df_count2$df.chr <- str_extract(df_count2$df.chr, "\\|.*\\|")
  df_count2$df.chr <- substr(df_count2$df.chr, 2, nchar(df_count2$df.chr) - 1)

  # merge with the dataset containing ref and alt of the DNPs
  data_info <- merge(df_count2, RefAltDNP, by.x = "df.chr", by.y = "POS")

  data_info2 <- data_info %>%
    mutate(is.alt = case_when(
      dfV5.2==alt ~ "ALT",
      dfV5.2!=alt ~ "NO",
    ))

  # create dataset with only reference genotypes
  ref_dt <- subset(data_info2,dfV5.2==".|." )
  colnames(ref_dt) <- c("DNPs","geno","n_ref","ref","alt","is.alt")

  # create dataset with only alternative genotypes
  alt_dt <- subset(data_info2,is.alt=="ALT")
  colnames(alt_dt) <- c("DNPs","geno","n_alt","ref","alt","is.alt")

  # create dataset with all other genotypes
  altro_dt <- subset(data_info2,is.alt!="ALT" & dfV5.2!=".|.")
  altro_dt2 <- data.frame(altro_dt$df.chr, altro_dt$n)

  altro_dt3 <- altro_dt2 %>%
    group_by(altro_dt.df.chr) %>%
    summarise(sum = sum(altro_dt.n))
  colnames(altro_dt3) <- c("DNPs","n_altro")

  # merge the datasets based on DNPs
  merged1 <- merge(ref_dt,alt_dt,by.x = "DNPs", by.y = "DNPs",all = T)
  merged2 <- merge(merged1,altro_dt3,by.x = "DNPs", by.y = "DNPs", all=T)
  merged3 <- merge(merged2,RefAltDNP,by.x = "DNPs", by.y = "POS", all=T)

  merged3[is.na(merged3)] <- 0

  dataframe_final <- merged3[,c(1,13,14,3,8,12)]
  colnames(dataframe_final) <- c("DNPs","Reference","Alternative","N_Reference", "N_Alternative", "N_Altro")

  # add column for REF / (REF + ALT)
  dataframe_final <- dataframe_final %>%
    mutate(ratio_Reference = N_Reference / (N_Reference + N_Alternative),
           total_reads = (N_Reference + N_Alternative + N_Altro))

  # filter out all DNPs where N_Altro > 10% of total_reads
  dataframe_final <- dataframe_final %>%
    filter(N_Altro / total_reads <= 0.10)

  filtered_2 <- dataframe_final %>%
    mutate(New_Column = ifelse(ratio_Reference == 1, ".|.", Alternative))

  filtered_3 <- merge(filtered_2, altro_dt, by.x = "DNPs", by.y = "df.chr")
  filtered_3 <- filtered_3[, -c(12:14)]

  # separate observed bases, reference bases, and alternative bases into two columns each
  filtered_4 <- filtered_3 %>%
    mutate(
      observed_base1 = str_split_fixed(dfV5.2, "\\|", 2)[, 1],   # Prima base osservata
      observed_base2 = str_split_fixed(dfV5.2, "\\|", 2)[, 2],   # Seconda base osservata

      ref_base1 = str_split_fixed(Reference, "\\|", 2)[, 1],     # Prima base di Reference
      ref_base2 = str_split_fixed(Reference, "\\|", 2)[, 2],     # Seconda base di Reference

      alt_base1 = str_split_fixed(Alternative, "\\|", 2)[, 1],   # Prima base di Alternative
      alt_base2 = str_split_fixed(Alternative, "\\|", 2)[, 2]    # Seconda base di Alternative
    )


  # create filtered_3a containing only rows where ratio_Reference equals 1
  filtered_3a <- filtered_4 %>%
    filter(ratio_Reference == 1)

  # create filtered_3b containing only rows where ratio_Reference equals 0
  filtered_3b <- filtered_4 %>%
    filter(ratio_Reference == 0)


  #### Errors for homozygous REF ####
  filtered_3a <- filtered_3a %>%
    mutate(
      single_error_b1 = ifelse(observed_base1 != "." & observed_base2 == ".", n, 0),
      single_error_b2 = ifelse(observed_base1 == "." & observed_base2 != ".", n, 0),
      double_error = ifelse(observed_base1 != "." & observed_base2 != ".", n, 0))

  summary_table_homRef <- filtered_3a %>%
    group_by(DNPs) %>%
    summarise(
      combinations = paste(unique(dfV5.2), collapse = "-"),
      total_single_error_b1 = sum(single_error_b1, na.rm = TRUE),
      total_single_error_b2 = sum(single_error_b2, na.rm = TRUE),
      total_double_error = sum(double_error, na.rm = TRUE))


  #### Errors for homozygous ALT ####
  filtered_3b <- filtered_3b %>%
    mutate(
      single_error_b1 = ifelse(observed_base1 != alt_base1 & observed_base2 == alt_base2, n, 0),
      single_error_b2 = ifelse(observed_base1 == alt_base1 & observed_base2 != alt_base2, n, 0),
      double_error = ifelse(observed_base1 != alt_base1 & observed_base2 != alt_base2, n, 0))

  summary_table_homAlt <- filtered_3b %>%
    group_by(DNPs) %>%
    summarise(
      combinations = paste(unique(dfV5.2), collapse = "-"),
      total_single_error_b1 = sum(single_error_b1, na.rm = TRUE),
      total_single_error_b2 = sum(single_error_b2, na.rm = TRUE),
      total_double_error = sum(double_error, na.rm = TRUE))

  filtered_def <- bind_rows(summary_table_homRef,summary_table_homAlt)
  df_def <- merge(dataframe_final, filtered_def, by = "DNPs")


  # ERROR RATE estimation
  df_definitivo <- df_def %>%
    mutate(
      error_rate1 = total_single_error_b1 / total_reads,
      error_rate2 = total_single_error_b2 / total_reads,
      double_error_rate = 9/8*(total_double_error / total_reads))

  df_definitivo <- df_definitivo %>%
    rename(
      N_Other = N_Altro
      )

  df_summary <- df_definitivo %>%
    summarise(
      mean_error_rate1 = mean(error_rate1, na.rm = TRUE),
      sd_error_rate1 = sd(error_rate1, na.rm = TRUE),

      mean_error_rate2 = mean(error_rate2, na.rm = TRUE),
      sd_error_rate2 = sd(error_rate2, na.rm = TRUE),

      mean_double_error_rate = mean(double_error_rate, na.rm = TRUE),
      sd_double_error_rate = sd(double_error_rate, na.rm = TRUE)
    )

  # write final table and summary
  output_file <- file.path(output_dir, paste0("ErrorRate_DNPs_", base_name, ".txt"))
  summary_output_file <- file.path(output_dir, paste0("Summary_ErrorRate_DNPs_", base_name, ".txt"))

  write.table(df_definitivo, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(df_summary, file = summary_output_file, sep = "\t", row.names = FALSE, quote = FALSE)

}

for (file in files) {
  process_pileup_file(file)
}
