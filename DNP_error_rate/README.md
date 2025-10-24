# **DNP Error Rate**

DNP Error Rate is a post–processing module designed to estimate the **sequencing error rate** on Double Nucleotide Polymorphisms (DNPs) using `samtools mpileup` and a reference DNP list.  
It assumes that the pileup is generated exactly as required for **DNPcall** and does not process BAM files itself.

---

## **Reference**

**Letizia Pistacchia, Francesco Ravasini, Elisa Bella, Eugenia D’Atanasio, Fulvio Cruciani, Beniamino Trombetta**  
*DNPcall: A New Pipeline for Accurate Double Nucleotide Polymorphism Calling* Bioinformatics Advances, 2025; vbaf209  
https://doi.org/10.1093/bioadv/vbaf209

---

## **Software required**

- **samtools** ≥ 1.18 (for pileup generation)  
- **R** ≥ 4.2.0  
- **R packages**:
  - `optparse`
  - `tidyverse`
  - `readxl`

> All these packages are required by the script and must be installed before running.

---

## **Installation**

No installation is required. Just download the repository and move into the directory:

```bash
git clone https://github.com/lpistacchia/NEWPAT.git
cd NEWPAT/DNP_error_rate
```

---

## **Input requirements**

This script requires:

1) **Pileup file** generated with `samtools mpileup` exactly as required by **DNPcall**, i.e.:  
   - includes **1 flanking base upstream and 1 downstream** for each DNP  
   - produced following DNPcall pileup specification

2) **DNP list file** listing the two positions forming each DNP  
   (one row per DNP: pos1, pos2, ref, alt)

The script computes:
- **single–base error rate** (for each base composing the DNP)
- **overall DNP-level error rate**

---

## **Usage**

Run on a **single pileup file**:

```bash
Rscript compute_error_rate.R \
  --pileup path/to/sample.pileup.txt \
  --dnp_list path/to/DNP_list.txt \
  --output path/to/output_directory
```

---

## **Arguments**

| **Flag**           | **Description** |
|--------------------|-----------------|
| **`--pileup`**     | Path to a single samtools pileup file |
| **`--dnp_list`**   | Path to the DNP list file |
| **`--output`**     | Output directory (default: current folder) |

## **Output**

- **`ErrorRate_DNPs_<sample>.txt`**  
  Per–DNP counts and error rate estimates

- **`Summary_ErrorRate_DNPs_<sample>.txt`**  
  Mean and standard deviation of error rates across all DNPs

---

## **Notes**

- The script does **not** generate pileup files  
- Pileup must follow **DNPcall specification** and include **flanking bases**  
- Only listed DNPs are evaluated  
- Indels and inconsistent genotypes are filtered internally

---

## **Contact**

For questions or requests, please contact:

Letizia Pistacchia 
<letizia.pistacchia@uniroma1.it>

