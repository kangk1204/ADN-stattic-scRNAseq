# ADN-stattic-scRNAseq

This repository contains the R scripts used to reproduce the single-cell
RNA-seq analysis of mouse kidney samples from three conditions:

- `sham`: control kidney
- `ADN`: adenine-induced nephropathy
- `ADN_Stattic`: adenine-induced nephropathy treated with Stattic, a STAT3 inhibitor

The workflow starts from Cell Ranger output folders, performs quality control
and integration, annotates major kidney cell populations, and generates the
figure panels used in the manuscript.

## What is included

```
scripts/
  01_QC_and_Integration.R                 # QC, SoupX, doublet detection, Harmony integration
  02_Figure3_FigureS3.R                   # Figure 3A-F and Figure S3B-D
  03_Figure4.R                            # Figure 4A-C
  04_FigureS8.R                           # Figure S8A-D
  05_Myeloid_Macro_Figure5_S9.R           # Figure 5A-D and Figure S9A-B
  06_CellChat_Figure5E_S9C.R              # Figure 5E and Figure S9C
  07_Fibroblast_subset_FigureS10ABC.R     # Figure S10A-C
  08_MacroFibro_CellChat_Figure5F_S10DE.R # Figure 5F and Figure S10D-E
```

Raw sequencing data, Cell Ranger output folders, intermediate `.rds` files, and
generated figures are not stored in this repository.

## Before you start

You need three Cell Ranger sample folders with this structure:

```
HN00216204_10X_RawData_Outs/
  sham/
    outs/
      raw_feature_bc_matrix/
      filtered_feature_bc_matrix/
  ADN/
    outs/
      raw_feature_bc_matrix/
      filtered_feature_bc_matrix/
  ADN_Stattic/
    outs/
      raw_feature_bc_matrix/
      filtered_feature_bc_matrix/
```

`SoupX::load10X()` reads both the raw and filtered matrices, so each full
`outs/` folder is required for script `01_QC_and_Integration.R`.

You also need the gene-set text files used for module scoring and heatmaps.
These files are exported from `Table_S1_Gene_List.xlsx`, one sheet per
tab-delimited `.txt` file, with one gene symbol per row.

Expected gene-set files:

```
LipidOxidation.txt
EpigeneticRegulation.txt
InflammatoryResponse.txt
ECMOrganization.txt
Ferroptosis.txt
group1.txt
group2.txt
group3.txt
group4.txt
final_PST_Heatmap.txt
group4_GeneList.txt
group4_InflammatoryResponse.txt
group4_ECM.txt
HALLMARK_IL6_JAK_STAT3_SIGNALING.txt
```

## Recommended folder layout

The scripts can run from any location as long as the paths are set correctly.
For a simple setup, keep your analysis files in one project directory:

```
ADN_stattic_analysis/
  cellranger/
    HN00216204_10X_RawData_Outs/
  results/
    RDS/
    Figures/
    GeneSets/
```

Place the exported gene-set `.txt` files in `results/GeneSets/`.

## R environment

The analysis was tested with:

```
R 4.4.2
Seurat 5.3.0
SoupX 1.6.2
scDblFinder 1.20.2
harmony 1.2.3
scCustomize 3.0.1
scProportionTest 0.0.0.9
EnhancedVolcano 1.22.0
CellChat 2.1.2
FSA 0.10.1
ggsignif 0.6.4
pheatmap 1.0.13
```

Install packages from CRAN, Bioconductor, or the package source used by your
analysis environment. Some packages, such as `CellChat`, `scCustomize`, and
`scProportionTest`, may require installation from GitHub depending on your R
setup.

Every script writes `sessionInfo()` to the RDS output directory so that package
versions can be checked after the run.

## How to run the workflow

Open each script and edit the configuration block near the top. These are the
paths you usually need to change:

```r
data_root     <- "path/to/cellranger/HN00216204_10X_RawData_Outs"  # script 01 only
out_dir       <- "path/to/results/RDS"
fig_dir       <- "path/to/results/Figures"
gene_list_dir <- "path/to/results/GeneSets"
```

Create the output folders before running the scripts, or let the scripts create
them when `dir.create()` is already included.

Run the scripts in numerical order:

```r
source("scripts/01_QC_and_Integration.R")
source("scripts/02_Figure3_FigureS3.R")
source("scripts/03_Figure4.R")
source("scripts/04_FigureS8.R")
source("scripts/05_Myeloid_Macro_Figure5_S9.R")
source("scripts/06_CellChat_Figure5E_S9C.R")
source("scripts/07_Fibroblast_subset_FigureS10ABC.R")
source("scripts/08_MacroFibro_CellChat_Figure5F_S10DE.R")
```

The scripts share the same `out_dir` and pass intermediate `.rds` files to the
next steps. If a later script cannot find an input file, check that the previous
script finished successfully and that the same `out_dir` is used throughout.

## Workflow map

```
01 -> Kidney_HarmonyIntegrated_dim30_res1.rds
       |
       +-> 02 -> Stattic_annotated.rds
                  |
                  +-> 03 -> PST_subset.rds  ----+
                  +-> 04 -> GEC_subset.rds  ----+--> 06
                  +-> 05 -> Myeloid_subset.rds  |
                  |        Macro_subset.rds  ---+--> 08
                  +-> 07 -> Fibro_subset.rds ---+
```

## Main outputs

- Integrated and annotated Seurat objects in `out_dir`
- Subset Seurat objects for PST, GEC, myeloid/macrophage, and fibroblast analyses
- CellChat objects for selected cell-cell communication analyses
- Figure files in `fig_dir`
- `sessionInfo_<script>.txt` files recording the R session used for each step

## Reproducibility notes

The scripts set `set.seed(1234)` before stochastic steps such as UMAP,
clustering, and doublet detection. Small differences can still occur if package
versions, BLAS/LAPACK libraries, or operating systems differ. If cluster numbers
shift between machines, use the marker panels in the annotation sections to
verify the intended cell identities before continuing.

## Troubleshooting

If `SoupX::load10X()` fails, confirm that `data_root` points to the directory
containing `sham`, `ADN`, and `ADN_Stattic`, and that each sample contains a full
Cell Ranger `outs/` folder.

If a module-score or heatmap step fails, check that all required gene-set `.txt`
files are present in `gene_list_dir` and that the first column contains gene
symbols.

If a later script reports that an `.rds` file is missing, rerun the earlier
script that creates it and confirm that all scripts use the same `out_dir`.

## Citation

If you use this code, please cite the associated manuscript when it becomes
available.
