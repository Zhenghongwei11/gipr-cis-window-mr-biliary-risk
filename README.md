# Drug-target MR: incretin pathway and biliary outcomes (reproducibility package)

This repository provides a lightweight reproduction bundle for a drug-target Mendelian randomization (MR) analysis of **incretin-pathway targets** (GLP1R, GIPR; comparator SLC5A2) and **biliary outcomes**, with pancreatitis retained as a secondary contextual endpoint.

## What you can reproduce quickly
Without any API tokens, you can regenerate the reported figures from the included small Source Data tables:

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
python3 -m pip install -r requirements.txt
```

## Notes
- The included `docs/source_data/` tables are traceable exports supporting the reported figures and tables.
- Full end-to-end recomputation from raw GWAS summary statistics is outside the scope of this lightweight bundle because it requires larger downloads and authenticated API access for some upstream queries.
