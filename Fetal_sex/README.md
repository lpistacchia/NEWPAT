# Fetal sex determination

This NEWPAT module determines fetal sex from maternal
cell-free DNA (cfDNA) using autosomal DNPs, X-DNPs
and Y-SNPs.

Paternal DNA is not required.

## Requirements

- R
- optparse
- dplyr

## Input

Three calling files from the same maternal cfDNA sample:

- Autosomal DNPcall output
- Chromosome X DNPcall output
- Chromosome Y SNPcall output

X and Y markers must be previously validated for specificity.

## Usage

```bash
Rscript Fetal_sex.R \
  --autosomes autosomal_DNPcall.txt \
  --X X_DNPcall.txt \
  --Y Y_SNPcall.txt \
  --out results
```

## Method

For autosomal and X-DNPs, the mean number of
non-maternal allele reads is calculated across
selected informative sites (MFF_A and MFF_X).

For Y-SNPs, the mean number of REF + ALT reads
is calculated across selected Y markers (MFF_Y).

Three ratios are calculated:

- A = MFF_X / MFF_A
- B = MFF_Y / MFF_A
- C = MFF_Y / MFF_X

Fetal sex is classified as F, M or Undetermined
using configurable thresholds.

## Parameters

| Parameter | Default | Description |
|---|---:|---|
| `--max_other_fraction` | 0.10 | Maximum Other read fraction |
| `--maternal_hom_threshold` | 0.90 | Maternal homozygosity threshold |
| `--min_fetal_reads` | 2* | Minimum non-maternal reads |
| `--min_fetal_percent` | Not set | Minimum non-maternal read percentage |
| `--include_zero_coverage_Y` | TRUE | Include Y markers with zero REF + ALT coverage |
| `--A_threshold` | 0.5 | A classification threshold |
| `--B_female_max` | 0.1 | Maximum B for female classification |
| `--B_male_min` | 1 | Minimum B for male classification |
| `--C_female_max` | 0.1 | Maximum C for female classification |
| `--C_male_min` | 2 | Minimum C for male classification |

*The default of 2 reads applies when neither fetal
read filter is specified. When both filters are
provided, both conditions must be satisfied.

## Classification

| Ratio | Female | Male |
|---|---:|---:|
| A | >= 0.5 | < 0.5 |
| B | <= 0.1 | >= 1 |
| C | <= 0.1 | >= 2 |

All three conditions must be satisfied for classification.
Otherwise, the result is Undetermined.

Classification thresholds are preliminary and require
further validation.

## Output

- `Fetal_sex_summary.txt`
- `Autosomal_informative_sites.txt`
- `X_informative_sites.txt`
- `Y_sites.txt`

## Examples

```bash
Rscript Fetal_sex.R \
  --autosomes examples/inputs/Female_fetus_autosomes.txt \
  --X examples/inputs/Female_fetus_X.txt \
  --Y examples/inputs/Female_fetus_Y.txt \
  --out example_results
```

## References

Pistacchia, L., et al. (2025).
DNPcall: A new pipeline for accurate double nucleotide
polymorphism calling. Bioinformatics Advances.
https://doi.org/10.1093/bioadv/vbaf209

D’Atanasio, E., et al. (2026).
A new noninvasive prenatal paternity test based on
double nucleotide polymorphisms abates sequencing
uncertainties and improves sensitivity. Genome Research.
https://doi.org/10.1101/gr.281702.125

## Contact

letizia.pistacchia@uniroma1.it
