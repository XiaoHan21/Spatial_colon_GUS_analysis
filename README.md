# Spatial_colon_GUS_analysis

Code for the spatial transcriptomic analysis of mouse colon tissue presented in **Extended Data Fig. 3** of the *Nature* article, entitled "*Gut Microbiota Mediate the Metabolism of
Colonic Prostaglandins*".

> **Extended Data Fig. 3 | Spatial transcriptomic analysis of mouse colon tissue revealed that genes encoding PGE2-Acyl-GlcA–processing GUS enzymes were enriched in the mucus-associated region.**

GF C57BL/6 mice were orally gavaged with a complex stool-derived microbial community prepared from SPF C57BL/6 donor mice. After 4 weeks, when the gut microbiota in the recipient mice had stabilized, colon tissues were collected and subjected to spatial transcriptomic analysis.

---

## Overview

| Panel | Content |
|:-----:|---------|
| **a** | Violin plots comparing spot-to-host distances between G1 (with a mean distance of ~130 μm from the host epithelial boundary, corresponding to the mucus-associated region) and G2 (with a mean distance of ~446 μm from the epithelial boundary, corresponding to the lumen-associated region) spatial groups across proximal, middle, and distal colon sections. Box plots show the median (centre line), 25th and 75th percentiles (box bounds), and whiskers extending to the minimum and maximum values. Sample sizes (G1 and G2, respectively) were 22,976 and 8,663 microbial spots for the proximal colon, 31,674 and 8,462 for the middle colon, and 19,281 and 9,514 for the distal colon. Within each colon segment, spot-to-host distances were compared between G1 and G2 using a two-sided Wilcoxon rank-sum test at the spot level. **** P < 0.0001 for each comparison. |
| **b** | Percentages of spots expressing PGE2-Acyl-GlcA–processing GUS genes in G1 and G2 regions of the proximal, middle and distal colon sections from one mouse. For each region (G1 or G2) of each section, the percentage was calculated as the number of GUS-positive spots (defined as spots with at least one uniquely mapped read assigned to any of the 11 active GUS genes) divided by the total number of microbial spots in that region. Numbers of GUS-positive spots/total microbial spots (G1 and G2, respectively) were as follows: proximal, 176/22,976 and 37/8,663; middle, 127/31,674 and 20/8,462; distal, 94/19,281 and 30/9,514. |
| **c** | Pie chart showing the relative proportions of PGE2-Acyl-GlcA–processing GUS-derived deduplicated transcript molecules, quantified by unique molecular identifiers (UMIs) in the mouse colon. |
| **d** | H&E staining and spatial visualization of distance categories (≤150 µm and >150 µm), total GUS-derived UMI counts (nUMIs; number of unique molecular identifiers), and species-specific GUS gene expression. Scale bar, 500 µm. |

---

## Workflow

```
 Kraken-classified reads (FASTQ) + GUS gene reference (FASTA)
                         │
                         ▼
        s01.bwa_align.bwaAdj.sh      BWA-MEM alignment → sorted BAM
                         │
                         ▼
        s02.generate_spatial_matrix.py   UMI-deduplicated spot × GUS gene matrix (CSV)
                         │
                         ▼
        s03.SpatialDataAnalysis_Visualisation.ipynb   Analysis & plotting → Extended Data Fig. 3a–d
```

## Repository structure

```
Spatial_colon_GUS_analysis/
├── s01.bwa_align.bwaAdj.sh                      # Step 1: alignment to GUS genes (BWA-MEM + SAMtools)
├── s02.generate_spatial_matrix.py               # Step 2: BAM → spatial expression matrix
├── s03.SpatialDataAnalysis_Visualisation.ipynb  # Step 3: analysis and visualisation (R)
└── README.md
```

---

## Requirements

| Software | Version |
|----------|---------|
| BWA-MEM | v0.7.19 |
| SAMtools | v1.23.1 |
| Python | v3.9.23 |
| R | v4.4.3 |

**Python packages** (Step 2): `pysam`, `pandas`

**R packages** (Step 3): `ggplot2`, `patchwork`, `ggpubr`, `Seurat`, `<add others>`

A conda environment (`spatial_gus`) was used for Steps 1–2, for example:

```bash
conda create -n spatial_gus python=3.9 bwa samtools pysam pandas -c bioconda -c conda-forge
conda activate spatial_gus
```

---

## Usage

### Step 1 · Sequence alignment — `s01.bwa_align.bwaAdj.sh`

Aligns Kraken-classified, single-end reads to the GUS gene reference with BWA-MEM, then converts, sorts and indexes the alignment with SAMtools.

| | |
|---|---|
| **Input** | `D04707D6.classified.fq` — Kraken-classified reads; read headers carry spatial coordinates (`Cx`, `Cy`), cell ID, UMI and Kraken taxid |
| | `GUS_dna.fna` — GUS gene reference (FASTA) |
| **Output** | `mapped_GUS.sorted.bam` and `mapped_GUS.sorted.bam.bai` |

Processing steps:
1. Replace spaces with underscores in FASTQ headers (prevents header truncation in the BAM file).
2. Index the reference (`bwa index`).
3. Align with BWA-MEM: `-k 13 -A 1 -B 2 -O 5 -T 15`, 30 threads.
4. Convert SAM → BAM, sort and index (`samtools`).

The script was written for a PBS cluster (`#PBS` header). **Edit the paths at the top of the script** (`OUTDIR`, `FQ_IN`, `FA_REF`) and the resource settings to match your environment, then run:

```bash
# on a PBS cluster
qsub s01.bwa_align.bwaAdj.sh

# or directly
bash s01.bwa_align.bwaAdj.sh
```

### Step 2 · Spatial expression matrix — `s02.generate_spatial_matrix.py`

Parses each aligned read, applies a MAPQ filter, and counts **unique UMIs** per GUS gene for each spatial spot (identified by `CellID`, `Cx`, `Cy`).

| | |
|---|---|
| **Input** | `mapped_GUS.sorted.bam` from Step 1 |
| **Output** | `<prefix>_MAPQ<threshold>.csv` |

Output format: one row per spot (`CellID`, `Cx`, `Cy`), one column per GUS gene, values are UMI counts after deduplication (0 if not detected).

| Argument | Description | Default |
|----------|-------------|---------|
| `-i`, `--input` | Input sorted BAM file (required) | — |
| `-o`, `--outprefix` | Output file prefix | `Spatial_GUS_Matrix` |
| `-q`, `--mapq` | Minimum MAPQ threshold | `10` |

The analysis in the manuscript used **MAPQ ≥ 15**:

```bash
python s02.generate_spatial_matrix.py \
    -i mapped_GUS.sorted.bam \
    -o Spatial_GUS_Matrix \
    -q 15
# → Spatial_GUS_Matrix_MAPQ15.csv
```

### Step 3 · Analysis and visualisation — `s03.SpatialDataAnalysis_Visualisation.ipynb`

An R notebook (R v4.4.3) that loads the matrices from Step 2 together with spot metadata, classifies spots into **G1** (mucus-associated) and **G2** (lumen-associated) groups based on their distance to the host epithelial boundary, and generates **Extended Data Fig. 3a–d**.

- **Input:** `<GUS matrices from Step 2>`, `Seurat objects>`
- **Output:** Extended Data Fig. 3a–d

Statistical tests: Wilcoxon rank-sum test (panel a) and Pearson's chi-squared test (panel b).

---

## Data availability

The supporting data are available on Zenodo ([DOI: 10.5281/zenodo.24333034](https://doi.org/10.5281/zenodo.24333034)). Further details are available from the contact below.

## Contact

`Xiao Han` · `han_xiao@ntu.edu.sg`
