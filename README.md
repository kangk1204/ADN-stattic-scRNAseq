# ADN-stattic-scRNAseq

Analysis code and gene lists accompanying the revised Advanced Science
manuscript on renal injury and Stattic treatment. The repository includes mouse
kidney single-cell RNA-seq, bulk RNA-seq gene set variation analysis (GSVA), KPMP
human kidney single-cell/single-nucleus RNA-seq, and CosMx spatial transcriptomics.

ADN denotes adenine-induced nephropathy, and Stattic is a STAT3 inhibitor.
Mouse single-cell conditions are `sham` (control), `ADN` and `ADN_Stattic`.

Expression datasets, processed Seurat objects, intermediate checkpoints and
generated figures are hosted separately or created during analysis. Experimental
images, Western blots and genome-browser tracks are outside this code collection.

## Analysis scripts

Run the numbered `.R` files as complete scripts. Figure references follow the
revised manuscript.

| Folder | Script | Analysis / panels |
|---|---|---|
| `scRNA-seq` | `01_QC_and_Integration.R` | SoupX, doublet detection, QC and Harmony integration; S3A |
| | `02_Figure3_FigureS3.R` | Atlas annotation, module scores and cell representation; 3A–F, S3B–D |
| | `03_Figure4_FigureS7.R` | PST transcriptional analysis; 4A–B, S7A; gene lists for S7B |
| | `04_FigureS8.R` | GEC subtypes, expression and cell representation; S8A–D; gene lists for S8E |
| | `05_Myeloid_Macro_Figure5_S10_S13.R` | Myeloid/macrophage subsets and expression programs; 5A–D, S10A–B, S13B |
| | `06_CellChat_Figure5E_S9_S10.R` | Renal-cell-to-myeloid communication; 5E, S9A–C, S10C–D |
| | `07_Fibroblast_subset_FigureS11ABC.R` | Fibroblast subsets, markers and cell representation; S11A–C |
| | `08_MacroFibro_CellChat_Figure5F_S11DEF.R` | Macrophage-to-fibroblast communication; 5F, S11D–F |
| `GSVA` | `01_Figure2E_GSVA.R` | Three gene-set scores across five experimental groups; 2E |
| `KPMP` | `01_KPMP_Atlas_FigureS4.R` | Cohort selection, atlas views and lineage/disease markers; S4 |
| | `02_KPMP_PT_S3_Figure4JK.R` | Human PT-S3 analysis; 4J–K |
| | `03_KPMP_GEC_Figure4LM.R` | Human GEC analysis; 4L–M |
| | `04_KPMP_Macrophage_FigureS12BC.R` | Human monocyte/macrophage analysis; S12B–C |
| `CosMx` | `01_CosMx_Import_and_SpecimenQC.R` | Import and audit nine biopsy specimens and spatial coordinates |
| | `02_Figure6C_Spatial_Maps.R` | Representative niche, cell-type and transcript-molecule maps; 6C |
| | `03_Figure6D_Specimenwise_Neighborhood.R` | Within-specimen distance and radial enrichment; 6D |

Gene lists are included in `scRNA-seq/gene_sets/` and `GSVA/gene_sets/`.
The mouse folder also contains `gene_sets.tsv`, `gene_set_manifest.tsv` and
`SHA256SUMS.tsv`; these define membership, aliases and checksums used by the
scripts. Keep them together. External pathway-enrichment results are needed for
panels whose enrichment calculation is outside the supplied code.

## Software requirements

The scripts specify the following reference versions:

| Software | Version |
|---|---|
| R | 4.4.2 |
| Seurat | 5.3.0 |
| GSVA | 1.48.3 |
| SoupX | 1.6.2 |
| scDblFinder | 1.20.2 |
| harmony | 1.2.3 |
| CellChat | 2.1.2 |
| scProportionTest | 0.0.0.9000 |

Additional dependencies include SeuratObject, Matrix, SingleCellExperiment,
BiocParallel, DelayedMatrixStats, future, ggplot2, patchwork, pheatmap, rstatix,
ggpubr, ggnewscale, digest, RANN, spatstat.geom and spatstat.explore. Each script
loads or checks the subset it needs. `ragg` is optional for additional CosMx PNG
exports. Install dependencies before running; scripts do not install packages.

Reference package sources:

- [CellChat source](https://github.com/jinworks/CellChat/tree/a0d3b2d231d46c8787177fffeac908270c253747)
- [scProportionTest source](https://github.com/rpolicastro/scProportionTest/tree/37a04900be1c991428da15af7f9aa7b0ad84661e)
- [GSVA 1.48.3 source archive](https://bioconductor.org/packages/3.17/bioc/src/contrib/GSVA_1.48.3.tar.gz)

GSVA 1.48.3 and its dependencies, including GSEABase, must be installed explicitly
if the current Bioconductor release selects another GSVA version. The bulk script
requires the 1.48.3 API. Several single-cell scripts also enforce their core
package versions. Check `packageVersion()` before execution; use an isolated R
library when maintaining the reference environment.

## Data and input files

Obtain input data separately. Access follows the study data-release schedule;
study accessions may remain private before manuscript acceptance.

| Workflow | Data source | Default local input |
|---|---|---|
| Bulk RNA-seq | [GSE329553](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329553) | `GSVA/data/kidney_GSVA.csv` |
| Mouse scRNA-seq | [GSE329554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329554) | `scRNA-seq/data/cellranger/` or `scRNA-seq/data/processed/` |
| KPMP | [KPMP atlas](https://atlas.kpmp.org/) | `KPMP/data/HumanKidney.rds` |
| CosMx | [GSE253439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253439) | `CosMx/data/GSE253439_seurat_object_cosmx.rds` |

The KPMP study access date was September 24, 2025. The CosMx source analysis is
available in [human_kidney_multiomics](https://github.com/mreck1/human_kidney_multiomics).
These local names are script interfaces, not a listing of deposit filenames.
Set paths when downloaded files use different names. Use the same data release
and processed objects as the study; later releases may differ.

**Bulk GSVA.** The TPM CSV requires unique mouse gene symbols in column 1 and
25 finite, nonnegative expression columns: five samples each for `Sham`,
`5D ADN`, `5D ADN_Stattic`, `14D ADN` and `14D ADN_Stattic`. Retain the full
expression background. Repeated group headers and supported replicate suffixes
are accepted. For unrelated sample IDs, define `config$sample_map` in the script
with `source_column`, `source_header`, `sample_id` and `group`, using the sample
metadata. Expression columns start at 2. Groups are not inferred from positions.

**Mouse scRNA-seq.** Raw processing begins with Cell Ranger outputs. The parent
folder must contain `sham/outs/`, `ADN/outs/` and `ADN_Stattic/outs/`, each with
`raw_feature_bc_matrix/`, `filtered_feature_bc_matrix/` and `analysis/` for SoupX.
FASTQ alignment and reference construction precede this workflow.

Alternatively, provide annotated Seurat checkpoints in `scRNA-seq/data/processed/`.
Explicit object paths take precedence; otherwise the original annotated name is
preferred over the generic name:

| Object | Original annotated name | Generic checkpoint name |
|---|---|---|
| Atlas | `Stattic_after_levels2.rds` | `Stattic_annotated.rds` |
| GEC | `GEC_subset_after_levels.rds` | `GEC_subset.rds` |
| Myeloid | `Myeloid.rename.rds` | `Myeloid_subset.rds` |
| Macrophage | `Macro.rename.rds` | `Macro_subset.rds` |
| Fibroblast | `Fibro_subset_after_levels.rds` | `Fibro_subset.rds` |
| Joint PCT/PST | `PT_subset_after_levels.rds` | `PT_CellChat_subset.rds` |

Objects need joined, normalized RNA layers, `condition` labels `sham`, `ADN` and
`ADN_Stattic`, and the applicable annotation: `CellType`, `GEC_CellType`,
`Myeloid_CellType`, `Macro_CellType`, `Fibro_CellType` or `PT_CellType`. Plotted
objects require stored UMAP. Subset cells and conditions are checked against the
parent. A processed atlas alone does not replace annotated subsets.
`Stattic_after_levels.rds` is a different atlas annotation; a generic
`Macro_subset.rds` from an original working folder may be unannotated.

**KPMP.** Supply the consortium-processed atlas with normalized RNA data, variable
features, an existing `umap` reduction and metadata `condition.long`, `cell_type`
and `subclass.full`. The embedded mapping defines broad lineages.
Healthy/Reference, AKI and CKD are retained; COVID-AKI and DKD are excluded.
`donor_id_column` defaults to `donor_id`; use supplied donor identifiers.
Figure S4E is skipped by default. To plot it, set `enrichment_tsv` to an external
TSV with `term`, `group` and `p_adjusted`; group names combine a lineage with
`__Healthy_Reference`, `__AKI` or `__CKD`. The script does not calculate enrichment.

**CosMx.** The object needs image containers `biopsy_1`, `biopsy_2` and `biopsy_3`,
metadata `Specimen_ID`, `FOV` and `Annotation.Lvl2`, cell centroids and SCT/data
expression. Figure 6C also needs `Niche` annotations, segmentation boundaries and
transcript-molecule coordinates. Embedded crop coordinates depend on the original
image orientation and data release. Image containers encompass multiple specimens.

## Configuration and execution

Edit the Configuration block, or set these environment variables to absolute
paths. Relative paths resolve from the folder containing the script. Annotation
maps and crop coordinates are embedded; no separate configuration download is
required.

| Variable | Input / setting |
|---|---|
| `KIDNEY_TPM_CSV` | Bulk TPM CSV |
| `KIDNEY_KPMP_RDS` | Original KPMP atlas |
| `KIDNEY_COSMX_RDS` | Processed CosMx object |
| `KIDNEY_CELLRANGER_DIR` | Parent of the three Cell Ranger sample folders |
| `KIDNEY_MOUSE_PROCESSED` | Directory of annotated mouse checkpoints |
| `KIDNEY_RESULTS_DIR` | Output directory for the selected workflow |

Mouse objects can also be set with `KIDNEY_MOUSE_ATLAS_RDS`,
`KIDNEY_MOUSE_GEC_RDS`, `KIDNEY_MOUSE_MYELOID_RDS`, `KIDNEY_MOUSE_MACROPHAGE_RDS`,
`KIDNEY_MOUSE_FIBROBLAST_RDS` and `KIDNEY_MOUSE_PT_RDS`, or the `processed_rds`
entries in each script.

Use a separate results directory for each workflow; keep its setting consistent
across dependent scripts. Without an override, each analysis folder uses its own
`results/`. Run each complete script in a fresh R session. For example, from the
repository root:

```sh
KIDNEY_TPM_CSV="/path/to/kidney_GSVA.csv" \
KIDNEY_RESULTS_DIR="/path/to/results/GSVA" \
Rscript --vanilla GSVA/01_Figure2E_GSVA.R

KIDNEY_CELLRANGER_DIR="/path/to/cellranger" \
KIDNEY_RESULTS_DIR="/path/to/results/scRNA-seq" \
Rscript --vanilla scRNA-seq/01_QC_and_Integration.R

KIDNEY_KPMP_RDS="/path/to/HumanKidney.rds" \
KIDNEY_RESULTS_DIR="/path/to/results/KPMP" \
Rscript --vanilla KPMP/01_KPMP_Atlas_FigureS4.R

KIDNEY_COSMX_RDS="/path/to/GSE253439_seurat_object_cosmx.rds" \
KIDNEY_RESULTS_DIR="/path/to/results/CosMx" \
Rscript --vanilla CosMx/01_CosMx_Import_and_SpecimenQC.R
```

The files also support `source("path/to/script.R")` in a fresh R session.

- **Mouse:** run 01–08 in order for raw processing. With annotated inputs, begin
  at 02 or the required downstream stage. Stages 03, 04, 05 and 07 use the atlas
  from 02. Stage 06 uses atlas, GEC and macrophage objects and a joint PCT/PST
  analysis; the PST-only object from 03 cannot replace this input. Stage 08 uses
  macrophage, fibroblast, joint PT and GEC objects from 05, 07, 06 and 04.
- **KPMP:** run 01, then 02–04. The latter stages independently read
  `01_Atlas/KPMP_filtered_atlas.rds` from the shared KPMP results directory.
- **CosMx:** run 01, then 02 and/or 03. Stage 03 uses the compact checkpoint
  from 01 and does not depend on 02's figures.
- **GSVA:** run the single script after preparing TPM input and group mapping.

New mouse clustering requires marker-based annotation review. The stage exports
marker tables, an unannotated checkpoint, assignments and a review token, then
stops before applying labels. Review the marker evidence and embedded map, enter
the corresponding token in `reviewed_cluster_maps` in that script, set
`restart_failed_stage=TRUE`, and rerun. Tokens depend on assignments and mapping;
matching cluster counts alone is insufficient. Supplied annotated checkpoints
retain annotations and embeddings and bypass this new-cluster review.

Completed stage directories are protected from overwriting. Use a new results
directory for another completed analysis. `restart_failed_stage=TRUE` archives
only an incomplete stage before retrying after its input or configuration issue
has been resolved. Downstream scripts require completed upstream checkpoints or
configured external processed objects.

## Outputs and interpretation

Outputs include Seurat/CellChat checkpoints, count and gene-matching tables,
statistical results and figures. Where implemented, configuration, `sessionInfo`,
loaded package versions and input/code/output SHA-256 manifests accompany results.
A `SUCCESS` file marks completion. Exact reproduction depends on the data release,
sample mapping, annotations, gene membership and software environment as well as
random seeds. Large single-cell objects and CellChat inference need sufficient
RAM and disk space for intermediate objects.

Mouse scRNA-seq has one pooled library from two mice per condition: three
libraries in total. Cell-level tests and recovered-cell proportions describe
those preparations and do not estimate between-mouse variability. Bulk GSVA
uses five samples per group with within-gene-set ANOVA/Tukey comparisons. KPMP
cell-level expression tests are exploratory; exported donor medians are
descriptive. CellChat probabilities describe inferred ligand–receptor communication.

CosMx distance tests use paired specimen medians and Holm adjustment across target
populations. Cell violins display underlying distributions; radial-enrichment
curves are descriptive. Default observation windows are estimated from FOV bounds.
Spatial proximity and inferred communication do not establish physical binding or
causal signaling. Refer to the article for experimental procedures, data access
conditions and analyses outside this repository.

Please cite the associated article when using this code.
