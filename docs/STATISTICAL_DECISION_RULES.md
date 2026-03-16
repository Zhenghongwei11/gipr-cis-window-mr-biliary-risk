# Statistical decision rules

This document records the predeclared statistical thresholds used for the reported analyses and figure outputs.

## Effect estimation
- Primary estimator: inverse-variance weighted (IVW) for multi-variant instruments.
- Single-variant instruments: Wald ratio.

## Multiplicity control
- False discovery rate (FDR) control used the Benjamini–Hochberg procedure.
- Results are interpreted with attention to both nominal p values and FDR-adjusted q values.

## Reporting
- Effect sizes are reported as odds ratios (OR) with 95% confidence intervals for binary outcomes.
- PheWAS screening used a stringent threshold of p < 1×10^-8 to characterise non-target association burden.
