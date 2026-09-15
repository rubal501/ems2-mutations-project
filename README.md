# TP53 ESM-2 embedding analysis — code

Code supporting the final report: point mutants of the human TP53 DNA-binding
domain are embedded with ESM-2, and WT-vs-mutant distances (global and local
k=0/3/5) are compared against physicochemical annotations, MAVE-DB functional
scores and BLOSUM62.

## Layout

```
code/
├── requirements.txt        Python dependencies (R packages listed below)
├── notebooks/              Analysis narratives (run top-to-bottom)
├── scripts/                Standalone script export of the HF pipeline
├── r_analysis/             Statistical tests on the embedding-distance table
├── data/                   Input table consumed by the R scripts
└── results/                Exported figures
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

- `violin_plots_k5.pdf` — k5 cosine-distance distributions by group.

## Environment

```bash
pip install -r requirements.txt
```

R packages: `dplyr`, `car`, `MASS`, `ggplot2`, `effsize`
(`install.packages(c("dplyr", "car", "MASS", "ggplot2", "effsize"))`).

## Notes

- Notebooks download the TP53 sequence (UniProt `P04637`) and MAVE-DB scores at
  runtime, so they need network access.
- Figures used by the report are stored in `../figures/`; notebook/script output
  paths are relative to each file's working directory.
