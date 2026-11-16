# Spatial_colon_GUS_analysis

Code for the spatial transcriptomic analysis of mouse colon tissue presented in **Extended Data Fig. 4** of the *Nature* manuscript.

> **Extended Data Fig. 4 | Spatial transcriptomic analysis of mouse colon tissue revealed that genes encoding PGE2-Acyl-GlcA–processing GUS enzymes were enriched in the mucus-associated region.**

Germ-free (GF) C57BL/6 mice were orally gavaged with a complex stool-derived microbial community prepared from SPF C57BL/6 donor mice. After 4 weeks, once the gut microbiota had stabilized, colon tissues were collected for spatial transcriptomic analysis (Stereo-seq). This repository contains the pipeline used to map microbial reads to GUS genes, build spot-level GUS expression matrices, and generate Extended Data Fig. 4a–d.

---

## Overview

| Panel | Content |
|:-----:|---------|
| **a** | Violin plots of spot-to-host distances between G1 (mucus-associated, mean ~130 µm from the epithelial boundary) and G2 (lumen-associated, mean ~446 µm) groups in proximal, middle and distal colon |
| **b** | Percentage of GUS-expressing spots in the mucus- vs. lumen-associated regions |
| **c** | Taxonomic composition of GUS-derived UMIs by species |
| **d** | H&E staining and spatial maps of distance categories (≤150 µm / >150 µm), total GUS UMI counts (nUMIs), and species-specific GUS gene expression |

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
        s03.SpatialDataAnalysis_Visualisation.ipynb   Analysis & plotting → Extended Data Fig. 4a–d
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

An R notebook (R v4.4.3) that loads the matrices from Step 2 together with tissue images and spot metadata, classifies spots into **G1** (mucus-associated) and **G2** (lumen-associated) groups based on their distance to the host epithelial boundary, and generates **Extended Data Fig. 4a–d**.

- **Input:** `<GUS matrices from Step 2>`, `<H&E images>`, `<spot-to-host distance / Seurat objects>`
- **Output:** Extended Data Fig. 4a–d (PDF)

Statistical tests: Wilcoxon rank-sum test (panel a) and Pearson's chi-squared test (panel b).

---

## Data availability

The supporting data are available from the contact below upon reasonable request.

## Citation

If you use this code, please cite:

> `<Authors>. <Title>. Nature (<year>). <DOI>`

## License

`<e.g. MIT License — see LICENSE for details>`

## Contact

`Xiao Han` · `han_xiao@ntu.edu.sg`
