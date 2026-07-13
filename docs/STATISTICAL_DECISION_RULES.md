# Statistical decision rules

This document records the predefined statistical thresholds used for the reported analyses and figure outputs.

## Effect estimation
- Primary estimator: inverse-variance weighted (IVW) for multi-variant instruments.
- Single-variant instruments: Wald ratio.

## Multiplicity control
- False discovery rate (FDR) control used the Benjamini–Hochberg procedure.
- Results are interpreted with attention to both nominal p values and FDR-adjusted q values.

## Reporting
- Effect sizes are reported as odds ratios (OR) with 95% confidence intervals for binary outcomes.
- Variant association profiles used retrieved association records at p < 1×10^-8 and p < 1×10^-5 to describe non-target association patterns.
