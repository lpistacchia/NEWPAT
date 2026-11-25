# Fetal sex determination

This NEWPAT module determines fetal sex from maternal
cell-free DNA (cfDNA) using autosomal DNPs, X-DNPs
and Y-SNPs.

## Requirements

- R
- optparse
- dplyr

## Input files

Three calling outputs from the same maternal cfDNA sample
are required:

| Argument | Input |
|---|---|
| --autosomes | Autosomal DNPcall output |
| --X | X chromosome DNPcall output |
| --Y | Y chromosome SNPcall output |
| --out | Output directory |

X and Y markers must be selected for chromosome specificity
before running the analysis.

Male controls, alleged father samples and population allele
frequencies are not required by this module.

## Method

### 1. Autosomal and X-DNP selection

Markers are selected according to:

- Maximum Other / total_count fraction: 0.10.
- Maternal homozygosity threshold: 0.90.
- Minimum non-maternal read support: 2 reads by default.

Maternal homozygous REF sites are identified when
REF / (REF + ALT) >= 0.90.

Maternal homozygous ALT sites are identified when
REF / (REF + ALT) <= 0.10.

At homozygous REF sites, ALT reads are considered
non-maternal reads.

At homozygous ALT sites, REF reads are considered
non-maternal reads.

Non-maternal read percentage is calculated as:

    100 * non-maternal reads / (REF + ALT)

The minimum number and percentage of non-maternal reads
can be specified independently or simultaneously.

When both thresholds are specified, both must be satisfied.
When neither is specified, the minimum is 2 reads.

### 2. Y-SNP selection

Y chromosome read support is calculated as:

    Y_reads = N_REF + N_ALT

By default, all Y-SNPs are included, including markers
with zero REF + ALT reads.

Zero-read markers can optionally be excluded.

If no Y-SNP has REF or ALT reads, MFF_Y is set to zero.

### 3. Mean fetal read support

    MFF_A = mean non-maternal reads across selected autosomal DNPs

    MFF_X = mean non-maternal reads across selected X-DNPs

    MFF_Y = mean REF + ALT reads across included Y-SNPs

### 4. Fetal sex ratios

    A = MFF_X / MFF_A

    B = MFF_Y / MFF_A

    C = MFF_Y / MFF_X

### 5. Fetal sex classification

Default classification thresholds:

| Ratio | Female (F) | Male (M) |
|---|---|---|
| A | >= 0.5 | < 0.5 |
| B | <= 0.1 | >= 1 |
| C | <= 0.1 | >= 2 |

All three conditions must support the same classification.

Discordant or undefined ratios result in an Undetermined
classification.

Classification thresholds are configurable and should be
validated for the specific sequencing panel and dataset.

## Optional parameters

| Parameter | Default | Description |
|---|---|---|
| --max_other_fraction | 0.10 | Maximum Other / total_count fraction |
| --maternal_hom_threshold | 0.90 | Maternal homozygosity threshold |
| --min_fetal_reads | 2* | Minimum non-maternal reads |
| --min_fetal_percent | NULL | Minimum non-maternal read percentage |
| --include_zero_coverage_Y | TRUE | Include Y markers with zero REF + ALT reads |
| --A_threshold | 0.5 | X/autosome classification threshold |
| --B_female_max | 0.1 | Maximum B for female classification |
| --B_male_min | 1 | Minimum B for male classification |
| --C_female_max | 0.1 | Maximum C for female classification |
| --C_male_min | 2 | Minimum C for male classification |

*The default of 2 reads applies when neither fetal read
threshold is specified.

## Usage

```bash
Rscript Fetal_sex.R \
  --autosomes autosomes.txt \
  --X X_DNPs.txt \
  --Y Y_SNPs.txt \
  --out results
```

To exclude Y markers with zero REF + ALT reads:

```bash
Rscript Fetal_sex.R \
  --autosomes autosomes.txt \
  --X X_DNPs.txt \
  --Y Y_SNPs.txt \
  --include_zero_coverage_Y FALSE \
  --out results_no_zero_Y
```

## Examples

Two synthetic examples are provided, each containing
20 autosomal DNPs, 8 X-DNPs and 10 Y-SNPs.

Run the following commands from the Fetal_sex directory.

### Female example

```bash
Rscript Fetal_sex.R \
  --autosomes examples/inputs/Female_autosomes.txt \
  --X examples/inputs/Female_X.txt \
  --Y examples/inputs/Female_Y.txt \
  --out Female_test
```

### Male example

```bash
Rscript Fetal_sex.R \
  --autosomes examples/inputs/Male_autosomes.txt \
  --X examples/inputs/Male_X.txt \
  --Y examples/inputs/Male_Y.txt \
  --out Male_test
```

### Expected results

| Variable | Female | Male |
|---|---:|---:|
| Autosomal informative sites | 16 | 16 |
| X informative sites | 6 | 4 |
| Y sites included | 10 | 10 |
| MFF_A | 20 | 20 |
| MFF_X | 20 | 2.5 |
| MFF_Y | 0 | 30 |
| Ratio A | 1 | 0.125 |
| Ratio B | 0 | 1.5 |
| Ratio C | 0 | 12 |
| Fetal sex | F | M |

These synthetic examples demonstrate software functionality
and do not constitute clinical validation.

## Output files

| File | Description |
|---|---|
| Fetal_sex_summary.txt | Summary of selected markers, read support, ratios, parameters and fetal sex |
| Autosomal_informative_sites.txt | Selected autosomal DNPs |
| X_informative_sites.txt | Selected X-DNPs |
| Y_sites.txt | Y-SNPs, read support and inclusion status |

## References

Pistacchia, L., et al. (2025).
DNPcall: A new pipeline for accurate double nucleotide
polymorphism calling.
*Bioinformatics Advances*, vbaf209.
https://doi.org/10.1093/bioadv/vbaf209

D'Atanasio, E., et al. (2026).
A new noninvasive prenatal paternity test based on double
nucleotide polymorphisms abates sequencing uncertainties
and improves sensitivity.
*Genome Research*. Advance online publication.
https://doi.org/10.1101/gr.281702.125

## Contact

letizia.pistacchia@uniroma1.it
