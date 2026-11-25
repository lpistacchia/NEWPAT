# Fetal sex determination

This NEWPAT module estimates fetal sex from maternal
cell-free DNA (cfDNA) using autosomal DNPs, X-chromosome
DNPs and Y-chromosome SNPs.

The methodology is described in the Genome Research
manuscript:

*A new non-invasive prenatal paternity test based on
double nucleotide polymorphisms abates sequencing
uncertainties and improves sensitivity.*

## Input

- Autosomal DNPcall output from maternal cfDNA.
- X-chromosome DNPcall output from maternal cfDNA.
- Y-chromosome SNPcall output from maternal cfDNA.
- Optional X-chromosome DNPcall outputs from male
  genomic DNA controls.

## Method

Autosomal and X-chromosome informative sites are
selected using the same criteria as the NEWPAT
prenatal CPI module:

- Maternal REF allele fraction >= 0.90 or <= 0.10.
- Other read fraction <= 0.10.
- At least 2 non-maternal reads by default.

MFF_A and MFF_X are calculated as the mean number
of non-maternal reads across informative autosomal
and X-chromosome sites, respectively.

MFF_Y is calculated as the mean number of REF + ALT
reads across all Y-chromosome markers, including
markers with zero coverage.

## Fetal sex determination

Three ratios are calculated:

| Ratio | Formula | Female | Male |
|-------|---------|--------|------|
| A | MFF_X / MFF_A | >= 0.5 | < 0.5 |
| B | MFF_Y / MFF_A | <= 0.1 | >= 1 |
| C | MFF_Y / MFF_X | <= 0.1 | >= 2 |

Fetal sex is assigned when all three ratios support
the same classification. Otherwise, the result is
reported as Undetermined.

These are configurable default classification
thresholds implemented in NEWPAT.

## Optional X-chromosome specificity control

X-DNP specificity can be evaluated using male
genomic DNA controls.

By default, a marker is excluded when the minor
allele fraction exceeds 0.01 in at least two male
controls with REF + ALT coverage >= 100.

Markers with insufficient data for specificity
assessment can also be excluded optionally.

## Usage

    Rscript Fetal_sex.R \
      --autosomes autosomal_DNPcall.txt \
      --X X_DNPcall.txt \
      --Y Y_SNPcall.txt \
      --out results

With X-chromosome specificity control:

    Rscript Fetal_sex.R \
      --autosomes autosomal_DNPcall.txt \
      --X X_DNPcall.txt \
      --Y Y_SNPcall.txt \
      --male_controls male_controls.txt \
      --out results

The male_controls.txt file contains one path to
an X-chromosome DNPcall output per line.

## Optional parameters

| Parameter | Default |
|-----------|---------|
| --min_fetal_reads | 2 if neither fetal-read filter is specified |
| --min_fetal_percent | Not applied |
| --male_controls | Not applied |
| --exclude_X | Not applied |
| --exclude_untested_X | FALSE |
| --min_male_depth | 100 |
| --max_minor_fraction | 0.01 |
| --min_males_failed | 2 |
| --A_threshold | 0.5 |
| --B_female_max | 0.1 |
| --B_male_min | 1 |
| --C_female_max | 0.1 |
| --C_male_min | 2 |

Both fetal-read filters can be applied simultaneously.

## Output

- Fetal_sex_summary.txt
- Autosomal_informative_sites.txt
- X_informative_sites.txt
- Y_sites.txt
- X_excluded_markers.txt
- X_specificity_QC.txt (when male controls are provided)

## Requirements

R packages: optparse, dplyr.
