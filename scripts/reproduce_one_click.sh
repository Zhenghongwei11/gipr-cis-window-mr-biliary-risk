#!/usr/bin/env bash
set -euo pipefail

mkdir -p plots/causal/drug_target_mr
mkdir -p plots/publication/pdf plots/publication/png

Rscript scripts/mr/plot_workflow_diagram.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_primary_forest.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_wald_diagnostics.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_phewas_contrast.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_regional_pleiotropy_locus.R --outdir plots/causal/drug_target_mr
python3 scripts/mr/build_variant_trait_pathway_map.py

cp -f plots/causal/drug_target_mr/workflow_diagram.pdf plots/publication/pdf/Figure1.pdf
cp -f plots/causal/drug_target_mr/workflow_diagram.png plots/publication/png/Figure1.png
cp -f plots/causal/drug_target_mr/regional_pleiotropy_locus.pdf plots/publication/pdf/Figure2.pdf
cp -f plots/causal/drug_target_mr/regional_pleiotropy_locus.png plots/publication/png/Figure2.png
cp -f plots/causal/drug_target_mr/forest_primary.pdf plots/publication/pdf/Figure3.pdf
cp -f plots/causal/drug_target_mr/forest_primary.png plots/publication/png/Figure3.png
cp -f plots/causal/drug_target_mr/wald_diagnostics.pdf plots/publication/pdf/Figure4.pdf
cp -f plots/causal/drug_target_mr/wald_diagnostics.png plots/publication/png/Figure4.png
cp -f plots/causal/drug_target_mr/variant_trait_pathway_heatmap.pdf plots/publication/pdf/Figure5.pdf
cp -f plots/causal/drug_target_mr/variant_trait_pathway_heatmap.png plots/publication/png/Figure5.png
cp -f plots/causal/drug_target_mr/phewas_contrast.pdf plots/publication/pdf/FigureS1.pdf
cp -f plots/causal/drug_target_mr/phewas_contrast.png plots/publication/png/FigureS1.png

python3 - <<'PY'
import hashlib
from pathlib import Path

root = Path(".")
rows = []
for base, asset_type, pattern in [
    (root / "docs" / "source_data", "source_data", "*.tsv"),
    (root / "results", "result_table", "*.tsv"),
    (root / "plots" / "publication" / "pdf", "figure_pdf", "*.pdf"),
    (root / "plots" / "publication" / "png", "figure_png", "*.png"),
]:
    for p in sorted(base.rglob(pattern)):
        if p.name == "DATA_MANIFEST.tsv":
            continue
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        rows.append((asset_type, p.as_posix(), h))

text = "asset_type\tpath\tsha256\n" + "".join(
    f"{t}\t{path}\t{sha}\n" for t, path, sha in rows
)
Path("docs/DATA_MANIFEST.tsv").write_text(text, encoding="utf-8")
Path("results/DATA_MANIFEST.tsv").write_text(text, encoding="utf-8")
PY

echo "OK: regenerated figures under plots/publication/."
