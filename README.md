# OG-CAF

THBS2+ matrix CAFs facilitate lung adenocarcinoma progression.

This repository holds the CAF subtype scoring code used in the study of cancer-associated fibroblasts (CAFs) during human lung adenocarcinoma (LUAD) progression. The study uses single-cell RNA sequencing, spatial transcriptomics, and in vivo and in vitro models, and identifies a THBS2-high matrix CAF (mCAF) subtype linked to poor outcome, recurrence and invasion, acting through a THBS2-SDC4-EMT axis.

## Citation

Ren Y, Ming R, Zuo A, Liu S, Ba Y, Zhang Y, Chen Y, Pan T, Luo P, Cheng Q, Deng J, Yue Y, Xu H, Weng S, Han X, Zhou D, Liu Z. Cancer-associated fibroblasts drive lung adenocarcinoma progression via THBS2-mediated epithelial-mesenchymal transition. *Oncogene* 2025;44(44):4284-4297. doi: [10.1038/s41388-025-03569-9](https://doi.org/10.1038/s41388-025-03569-9)

The published article is titled differently from this repository's description; the correspondence rests on its abstract (THBS2-high mCAFs in LUAD progression).

## Contents

| File | Purpose |
| --- | --- |
| `Analysis.R` | Scores each cell for iCAF, mCAF and apCAF gene sets with Seurat `AddModuleScore` and draws violin plots of the three scores per CAF group |
| `CAF_markers.rda` | `CAF_markers`, a list of 23 marker gene sets collected from published CAF studies (for example `iCAF_Cell`, `myCAF_Cell`, `apCAF_Cell`, `Myofibro`, `iCAF_MC`, `mCAF_MC`) |
| `Reference.xlsx` | Source table for those gene sets: year, journal, PMID, CAF type and genes |

## How the scores are built

`Analysis.R` forms three gene sets by merging the published markers with a short list of canonical genes written in the script:

- mCAF: `myCAF_Cell`, `Myofibro`, `mCAF_MC` and the canonical mCAF genes
- iCAF: `iCAF_Cell`, `iCAF_MC` and the canonical iCAF genes
- apCAF: `apCAF_Cell` and the canonical apCAF genes

It then calls `AddModuleScore(sc, features = list(mCAF, iCAF, apCAF))` and plots the scores with `VlnPlot` and `patchwork`.

## Running it

The script expects a Seurat object named `sc` of CAFs that carries a `CAF.type` column in `sc@meta.data`; the single-cell data and this object are not included in the repository. The score columns are renamed by position (`meta.data[37:39]`), so adjust those indices if your object has a different number of metadata columns. Run from the repository directory so that `CAF_markers.rda` is found. Required packages: Seurat (the script attaches `SCOPE`, then calls Seurat functions), ggplot2 and patchwork.

## License

GPL-3.0, see `LICENSE`.
