# SNPcall

SNPcall calls biallelic single nucleotide polymorphisms (SNPs) from BAM files.
It adapts the read-based workflow described for
[DNPcall](https://doi.org/10.1093/bioadv/vbaf209) to a single variant position
per marker.

## Requirements

- samtools ≥ 1.18
- R ≥ 4.2.0
- R packages `optparse`, `tidyverse`, and `ComplexUpset`

## Installation

```bash
git clone https://github.com/lpistacchia/NEWPAT.git
cd NEWPAT/SNPcall
```

## Input files

- **BAM list.** A text file with one BAM path per line. Paths can be absolute
  or relative to the directory containing the BAM list.
- **SNP list.** A tab-separated file with the header `chr`, `pos`, `ref`, `alt`.
  Positions are 1-based. REF and ALT must be distinct single nucleotides.
- **Indexed reference FASTA.** Use the same genome build and chromosome names
  as in the BAM files and SNP list.

The [example SNP list](examples/inputs/SNPs_list.txt) contains 15 SNPs.
The first two entries are shown below.

```tsv
chr	pos	ref	alt
chr1	100687460	C	A
chr1	103400491	G	A
```

## Usage

```bash
Rscript SNPcall_latest.R \
  --bamlist=examples/inputs/bamlist.txt \
  --SNPs=examples/inputs/SNPs_list.txt \
  --reference=path/to/reference.fa \
  --out=results_SNP
```

Replace `sample.bam` in [`examples/inputs/bamlist.txt`](examples/inputs/bamlist.txt)
with the full path to your BAM file. Use one BAM path per line. Replace
`path/to/reference.fa` with the path to your indexed reference FASTA.

| Argument | Description |
|---|---|
| `--bamlist` | Path to the BAM list |
| `--SNPs` | Path to the SNP list |
| `--reference` | Path to the indexed reference FASTA |
| `--out` | Path to a new output directory |
| `--parallel` | Process samples in parallel |
| `--cores` | Number of workers to use with `--parallel` |

Add `--parallel --cores=2` to process samples with two workers.

## Outputs

| Path | Contents |
|---|---|
| `individual_outputs/full_output_<sample>.txt` | Per-sample results, one row per SNP |
| `summary/<output_directory>_AllGenotypes.txt` | Genotype calls across samples |
| `summary/Unassigned_genotypes.pdf` | Number of unassigned genotypes per SNP |
| `summary/SNPs_covered_UpSet_plot.pdf` | SNPs with assigned genotypes across samples. Generated for runs with at least two samples. |
| `summary/resources.tsv` | Per-sample runtime and memory measurements |
| `plots/` | Per-sample read-count and genotype plots |
| `logs/` | Pileup generation and sample processing logs |

### Example outputs

The [example output table](examples/outputs/output_snpcall.txt) contains
the results for the 15 SNPs in the example panel.

- [Read-count and genotype plots](examples/outputs/plots/example_sample)
- [Unassigned genotypes plot](examples/outputs/summary/Unassigned_genotypes.pdf)

### Per-sample output columns

| Columns | Description |
|---|---|
| `chr`, `pos` | Chromosome and 1-based SNP position |
| `REF`, `ALT` | Reference and alternative nucleotides |
| `N_REF`, `N_ALT`, `N_Other` | Counts of reads matching REF, ALT, or another base |
| `total_count` | Sum of `N_REF`, `N_ALT`, and `N_Other` |
| `ratio` | REF fraction among REF and ALT reads, calculated as `N_REF / (N_REF + N_ALT)` |
| `GL0`, `GL1`, `GL2` | Natural-log genotype likelihoods for `0/0`, `0/1`, and `1/1`, respectively |
| `GT` | Genotype call (`0/0`, `0/1`, `1/1`, or `NA`) |
| `combined_bases`, `combined_counts` | Observed bases and corresponding read counts, listed in the same order |
| `reads_discarded_for_indels`, `ratio_reads_discarded_for_indels` | Count and proportion of read groups with an indel at the SNP position, calculated before indel filtering |

## Contact

For issues or questions about SNPcall, contact Letizia Pistacchia at
[letizia.pistacchia@uniroma1.it](mailto:letizia.pistacchia@uniroma1.it).

## Citation

Please cite the following publications if you use the NEWPAT pipeline:

- Pistacchia, L., et al. (2025). DNPcall: A new pipeline for accurate double
  nucleotide polymorphism calling. *Bioinformatics Advances*, vbaf209.
  https://doi.org/10.1093/bioadv/vbaf209
- D’Atanasio, E., et al. (2026). A new noninvasive prenatal paternity test
  based on double nucleotide polymorphisms abates sequencing uncertainties
  and improves sensitivity. *Genome Research*. Advance online publication.
  https://doi.org/10.1101/gr.281702.125
