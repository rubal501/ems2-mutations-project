# TP53 ESM-2 embedding analysis — code

## Introduction

Estimating the functional effect of missense mutations from sequence alone
remains a hard problem. Classical approaches such as BLOSUM62 capture aggregate
evolutionary evidence about amino-acid substitutions but ignore the sequence
context in which a mutation occurs, whereas protein language models (PLMs) like
ESM-2 produce contextual per-residue embeddings. This project asks whether a
*local* embedding representation centered on the mutated site captures missense
effects better than a *global* representation of the whole protein, and how
ESM-2 compares against a classic BLOSUM62 baseline.

TP53 — the most frequently mutated gene in human cancer — is an ideal testbed:
its DNA-binding domain (DBD) concentrates well-characterized mutational
hotspots, and large-scale MAVE functional scores are available for thousands of
variants (Giacomelli et al., MaveDB).

The analysis runs in two phases: (1) a synthetic panel of 228 point mutants
(12 DBD positions × 19 substitutions) is embedded with ESM-2, and WT–mutant
distances are computed globally and over local windows (k=0/3/5); (2) MAVE-DB
functional scores are integrated to test whether embedding-derived metrics
associate with measured function, and to benchmark ESM-2 against BLOSUM62 via
Spearman correlations and bootstrap confidence intervals.

Key findings: global averaging dilutes the effect of single substitutions, while
local (k=5) distances recover a stronger, more interpretable signal that
correlates significantly with experimental function. ESM-2 and BLOSUM62 show no
uniform winner: ESM-2 outperforms BLOSUM62 in the comparative DBD subset, while
BLOSUM62 remains competitive at canonical hotspots.

See [`manuscript/`](manuscript/) for the full report; this repository holds the
code, data and results behind it.

## Layout

```
repository/
├── manuscript/             Typst report sources (reporte.typ, anexo.typ, references.bib)
│   └── build/              Compiled PDFs (gitignored)
├── notes/                  Working notes (notas.md)
├── figures/                Figures used by the report
├── results/
│   ├── tables/             Generated CSV tables (pipeline outputs)
│   └── figures/            Exported figures
├── requirements.txt        Python dependencies (R packages listed below)
├── notebooks/              Analysis narratives (run top-to-bottom)
├── scripts/                Standalone script export of the HF pipeline
├── r_analysis/             Statistical tests on the embedding-distance table
└── data/                   Input table consumed by the R scripts
```

## Notebooks (`notebooks/`)

| File | Description |
| --- | --- |
| `esm2_pipeline_hf.ipynb` | **Canonical pipeline.** ESM-2 (`facebook/esm2_t33_650M_UR50D`) via Hugging Face `transformers`. Phases F1–F7: mutants, embeddings, distances, PCA/t-SNE/UMAP, physicochemical annotation + Mann–Whitney, MAVE-DB scores, ESM-2 vs BLOSUM62. |
| `esm2_pipeline_fair_esm.ipynb` | Same analysis using the `fair-esm` package (`esm2_t33_650M_UR50D`). Executed, with outputs. |
| `esm2_pipeline_fair_esm_src.ipynb` | Source-identical twin of the previous notebook, outputs cleared. |
| `esm2_extended_stats.ipynb` | Extended statistics: physicochemical flag tests, MAVE-DB correlations, BLOSUM62 benchmark, subgroup correlations and bootstrap CIs. |

> `esm2_pipeline_fair_esm.ipynb` and `esm2_pipeline_fair_esm_src.ipynb` share the
> same source; only the saved outputs differ.

## Scripts (`scripts/`)

- `embeddings_hf.py` — script export of `esm2_pipeline_hf.ipynb`. Writes its
  CSV/SVG outputs to the current working directory; run it from `scripts/`.

## R analysis (`r_analysis/`)

Both scripts read `../data/tp53_synthetic_mutants_embedding_metrics_annotated.csv`.

- `embedding_shift_analysis.R` — normality checks, ANOVA/Kruskal–Wallis,
  pairwise Wilcoxon (Holm), Cliff's delta by mutation group.
- `shift_chemical_groups.R` — Mann–Whitney tests of k5 cosine distance by
  physicochemical flag (charge, class, size, proline/glycine/cysteine, …).

## Data (`data/`)

- `tp53_synthetic_mutants_embedding_metrics_annotated.csv` — 228 single
  substitutions (12 DBD positions × 19 aa) with cosine/euclidean distances at
  global/k0/k3/k5 and physicochemical annotations. Input to the R scripts.

## Results (`results/`)

- `tables/` — generated CSV tables exported by the notebooks/scripts
  (`metric_comparison.csv`, `compare_rank.csv`, `mannwhitney_k5_results.csv`,
  `correlation-embedding-table.csv`, `esm_vs_BLOSUMboots_group.csv`, …).
- `figures/violin_plots_k5.pdf` — k5 cosine-distance distributions by group.

## Manuscript (`manuscript/`)

- `reporte.typ` — main report; `anexo.typ` — appendix; `references.bib` —
  bibliography. Both `.typ` files reference figures via `../figures/`.
  Compile with `typst compile --root .. reporte.typ build/reporte.pdf` from `manuscript/`.

## Environment

```bash
pip install -r requirements.txt
```

R packages: `dplyr`, `car`, `MASS`, `ggplot2`, `effsize`
(`install.packages(c("dplyr", "car", "MASS", "ggplot2", "effsize"))`).

## Notes

- Notebooks download the TP53 sequence (UniProt `P04637`) and MAVE-DB scores at
  runtime, so they need network access.
- Figures used by the report are stored in `figures/`; notebook/script output
  paths are relative to each file's working directory.
