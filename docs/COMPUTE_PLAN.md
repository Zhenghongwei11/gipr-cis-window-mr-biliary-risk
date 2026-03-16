# Compute plan

## Quick reproduction (default)
The included `docs/source_data/` tables are sufficient to regenerate the final manuscript figures on a typical laptop in minutes:

```bash
bash scripts/reproduce_one_click.sh
```

This quick path is designed to be deterministic and does not require any API credentials.

## Full recomputation (optional; not required for review)
Full recomputation from raw or API-served GWAS summary statistics can require large downloads and/or authenticated endpoints. It is intentionally not included in this lightweight bundle.
