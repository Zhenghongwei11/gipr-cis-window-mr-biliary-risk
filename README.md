# Drug-target MR: incretin pathway → biliary/pancreatic outcomes (reproducibility package)

This repository provides a lightweight, reviewer-facing reproduction bundle for a drug-target Mendelian randomization (MR) analysis of **incretin-pathway targets** (GLP1R, GIPR; comparator SLC5A2) and **biliary/pancreatic outcomes**.

## What you can reproduce quickly
Without any API tokens, you can regenerate the manuscript figures from the included small “Source Data” tables:

```bash
bash scripts/reproduce_one_click.sh
```

Outputs:
- `plots/publication/pdf/` (vector PDFs)
- `plots/publication/png/` (300 dpi PNGs)

## Environment
This repo includes `renv.lock` to pin R package versions.

```bash
Rscript -e 'install.packages("renv"); renv::restore(lockfile = "renv.lock")'
```

## Notes
- The included `docs/source_data/` tables are submission-facing, traceable exports supporting the reported figures.
- Full end-to-end recomputation from raw GWAS summary statistics is intentionally out of scope for this lightweight bundle.
