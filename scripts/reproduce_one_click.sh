#!/usr/bin/env bash
set -euo pipefail

mkdir -p plots/causal/drug_target_mr
mkdir -p plots/publication/pdf plots/publication/png

Rscript scripts/mr/plot_workflow_diagram.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_primary_forest.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_wald_diagnostics.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_phewas_contrast.R --outdir plots/causal/drug_target_mr
Rscript scripts/mr/plot_regional_pleiotropy_locus.R --outdir plots/causal/drug_target_mr

cp -f plots/causal/drug_target_mr/workflow_diagram.pdf plots/publication/pdf/Figure1.pdf
cp -f plots/causal/drug_target_mr/workflow_diagram.png plots/publication/png/figure1.png
cp -f plots/causal/drug_target_mr/forest_primary.pdf plots/publication/pdf/Figure2.pdf
cp -f plots/causal/drug_target_mr/forest_primary.png plots/publication/png/Figure2.png
cp -f plots/causal/drug_target_mr/wald_diagnostics.pdf plots/publication/pdf/Figure3.pdf
cp -f plots/causal/drug_target_mr/wald_diagnostics.png plots/publication/png/Figure3.png
cp -f plots/causal/drug_target_mr/regional_pleiotropy_locus.pdf plots/publication/pdf/Figure4.pdf
cp -f plots/causal/drug_target_mr/regional_pleiotropy_locus.png plots/publication/png/Figure4.png
cp -f plots/causal/drug_target_mr/phewas_contrast.pdf plots/publication/pdf/FigureS1.pdf
cp -f plots/causal/drug_target_mr/phewas_contrast.png plots/publication/png/FigureS1.png

echo "OK: regenerated figures under plots/publication/."
