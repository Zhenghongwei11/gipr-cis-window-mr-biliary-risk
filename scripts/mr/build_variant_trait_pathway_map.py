#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


OUT_TABLE = Path("docs/source_data/Source_Data_14_variant_trait_pathway_map.tsv")
OUT_PNG = Path("plots/causal/drug_target_mr/variant_trait_pathway_heatmap.png")
OUT_PDF = Path("plots/causal/drug_target_mr/variant_trait_pathway_heatmap.pdf")

PHEWAS = Path("docs/source_data/Source_Data_9_phewas_hits_rs17561351_rs10407429_p1e8.tsv")
INSTRUMENTS = Path("docs/source_data/Source_Data_5_GIPR_instruments_wide1mb.tsv")
WALD_DIR = Path("results/causal/drug_target_mr_runs/hba1c_wide1mb/diagnostics")

VARIANTS = {
    "rs10407429": "Retained in strict GIPR window",
    "rs17561351": "Excluded by strict GIPR window",
}

PHEWAS_TRAITS = [
    ("Glycaemic", "HbA1c instrument GWAS", "Glycated haemoglobin HbA1c levels (UKB data field 30750)"),
    ("Glycaemic", "Type 2 diabetes", "Type 2 diabetes"),
    ("Glycaemic", "Glucose", "Glucose levels (UKB data field 30740)"),
    ("Lipid/lipoprotein", "LDL cholesterol", "Low density lipoprotein cholesterol levels"),
    ("Lipid/lipoprotein", "HDL cholesterol", "High density lipoprotein cholesterol levels (UKB data field 30760)"),
    ("Lipid/lipoprotein", "Apolipoprotein B", "apolipoprotein B"),
    ("Lipid/lipoprotein", "Triglycerides", "triglycerides"),
    ("Hematologic", "Mean platelet volume", "Mean platelet volume"),
]

WALD_TRAITS = [
    ("Biliary", "Cholelithiasis", "wald_GIPR_cholelithiasis_primary.tsv"),
    ("Biliary", "Cholecystectomy FinnGen", "wald_GIPR_cholecystectomy_secondary.tsv"),
    ("Biliary", "Cholecystectomy UKB", "wald_GIPR_cholecystectomy_ukb.tsv"),
    ("Biliary", "Cholecystitis", "wald_GIPR_cholecystitis_finngen.tsv"),
    ("Pancreatic", "Acute pancreatitis", "wald_GIPR_acute_pancreatitis_finngen.tsv"),
]

EXPLICIT_MISSING = [
    ("Anthropometric", "BMI / adiposity", "No cached p<1e-8 hit for these exact HbA1c GIPR variants in Source Data 9"),
]


def direction(beta: float | None) -> str:
    if beta is None or math.isnan(beta):
        return "NA"
    if beta > 0:
        return "positive"
    if beta < 0:
        return "negative"
    return "zero"


def signed_logp(beta: float | None, pvalue: float | None) -> float | None:
    if beta is None or pvalue is None or math.isnan(beta) or math.isnan(pvalue):
        return None
    if pvalue <= 0:
        val = 300.0
    else:
        val = -math.log10(pvalue)
    return val if beta > 0 else -val if beta < 0 else 0.0


def best_phewas_row(phewas: pd.DataFrame, variant: str, trait_query: str) -> pd.Series | None:
    rows = phewas[
        (phewas["rsid"] == variant)
        & (phewas["trait"].str.lower() == trait_query.lower())
    ].copy()
    if rows.empty:
        # fall back to case-insensitive contains for traits with source-specific suffixes
        rows = phewas[
            (phewas["rsid"] == variant)
            & (phewas["trait"].str.lower().str.contains(trait_query.lower(), regex=False))
        ].copy()
    if rows.empty:
        return None
    rows["p_num"] = pd.to_numeric(rows["p"], errors="coerce")
    return rows.sort_values("p_num").iloc[0]


def add_row(rows: list[dict[str, object]], *, variant: str, domain: str, trait_label: str,
            source: str, beta: float | None, se: float | None, pvalue: float | None,
            evidence_type: str, notes: str) -> None:
    rows.append(
        {
            "variant": variant,
            "variant_status": VARIANTS[variant],
            "domain": domain,
            "trait_label": trait_label,
            "source": source,
            "evidence_type": evidence_type,
            "beta": "" if beta is None or math.isnan(beta) else f"{beta:.6g}",
            "se": "" if se is None or math.isnan(se) else f"{se:.6g}",
            "pvalue": "" if pvalue is None or math.isnan(pvalue) else f"{pvalue:.6g}",
            "direction": direction(beta),
            "signed_neg_log10_p": "" if signed_logp(beta, pvalue) is None else f"{signed_logp(beta, pvalue):.6g}",
            "notes": notes,
        }
    )


def build_table() -> pd.DataFrame:
    phewas = pd.read_csv(PHEWAS, sep="\t")
    instruments = pd.read_csv(INSTRUMENTS, sep="\t")
    rows: list[dict[str, object]] = []

    for variant in VARIANTS:
        inst = instruments[instruments["snp"] == variant].iloc[0]
        add_row(
            rows,
            variant=variant,
            domain="Glycaemic",
            trait_label="HbA1c instrument effect",
            source=str(INSTRUMENTS),
            evidence_type="instrument",
            beta=float(inst["beta"]),
            se=float(inst["se"]),
            pvalue=float(inst["pvalue"]),
            notes="HbA1c cis instrument association used for the GIPR broad-window proxy",
        )

        for domain, label, query in PHEWAS_TRAITS:
            # HbA1c instrument effect already came from the instrument table.
            if label == "HbA1c instrument GWAS":
                continue
            hit = best_phewas_row(phewas, variant, query)
            if hit is None:
                add_row(
                    rows,
                    variant=variant,
                    domain=domain,
                    trait_label=label,
                    source=str(PHEWAS),
                    evidence_type="phewas_p1e8_absent",
                    beta=None,
                    se=None,
                    pvalue=None,
                    notes=f"No cached p<1e-8 PheWAS hit matching '{query}'",
                )
                continue
            add_row(
                rows,
                variant=variant,
                domain=domain,
                trait_label=label,
                source=f"{PHEWAS}:{hit['id']}",
                evidence_type="phewas_p1e8",
                beta=float(hit["beta"]),
                se=float(hit["se"]) if pd.notna(hit["se"]) else None,
                pvalue=float(hit["p"]),
                notes=str(hit["trait"]),
            )

        for domain, label, note in EXPLICIT_MISSING:
            add_row(
                rows,
                variant=variant,
                domain=domain,
                trait_label=label,
                source=str(PHEWAS),
                evidence_type="not_observed_in_cached_p1e8_phewas",
                beta=None,
                se=None,
                pvalue=None,
                notes=note,
            )

        for domain, label, filename in WALD_TRAITS:
            path = WALD_DIR / filename
            wald = pd.read_csv(path, sep="\t")
            hit = wald[wald["snp"] == variant]
            if hit.empty:
                add_row(
                    rows,
                    variant=variant,
                    domain=domain,
                    trait_label=label,
                    source=str(path),
                    evidence_type="wald_absent",
                    beta=None,
                    se=None,
                    pvalue=None,
                    notes="Variant not present in this Wald diagnostic table",
                )
                continue
            rec = hit.iloc[0]
            add_row(
                rows,
                variant=variant,
                domain=domain,
                trait_label=label,
                source=str(path),
                evidence_type="wald_ratio",
                beta=float(rec["wald_beta"]),
                se=float(rec["wald_se"]),
                pvalue=float(rec["wald_pvalue"]),
                notes=f"Wald OR {float(rec['wald_or']):.3g} ({float(rec['wald_or_ci_lower']):.3g}-{float(rec['wald_or_ci_upper']):.3g})",
            )

    df = pd.DataFrame(rows)
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_TABLE, sep="\t", index=False)
    return df


def plot_heatmap(df: pd.DataFrame) -> None:
    plot_df = df.copy()
    plot_df["signed_neg_log10_p_num"] = pd.to_numeric(
        plot_df["signed_neg_log10_p"], errors="coerce"
    )
    order = [
        "HbA1c instrument effect",
        "Type 2 diabetes",
        "Glucose",
        "LDL cholesterol",
        "HDL cholesterol",
        "Apolipoprotein B",
        "Triglycerides",
        "Mean platelet volume",
        "BMI / adiposity",
        "Cholelithiasis",
        "Cholecystectomy FinnGen",
        "Cholecystectomy UKB",
        "Cholecystitis",
        "Acute pancreatitis",
    ]
    matrix = (
        plot_df.pivot(index="trait_label", columns="variant", values="signed_neg_log10_p_num")
        .reindex(order)
        .reindex(columns=list(VARIANTS))
    )

    values = matrix.to_numpy(dtype=float)
    finite = values[~pd.isna(values)]
    vmax = max(10.0, min(60.0, float(max(abs(finite))) if finite.size else 10.0))

    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad("#f2f2f2")
    im = ax.imshow(values, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(["rs10407429\n(strict retained)", "rs17561351\n(strict excluded)"], fontsize=9)
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8)
    ax.tick_params(length=0)

    for i, trait in enumerate(matrix.index):
        for j, variant in enumerate(matrix.columns):
            val = matrix.loc[trait, variant]
            label = "NA" if pd.isna(val) else f"{val:.1f}"
            color = "#555555" if pd.isna(val) else "white" if abs(val) > vmax * 0.45 else "#222222"
            ax.text(j, i, label, ha="center", va="center", fontsize=7, color=color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("signed -log10(P)", fontsize=9)
    fig.text(
        0.02,
        0.01,
        "Positive values indicate positive beta/Wald log-OR; negative values indicate inverse direction. Grey cells: no cached p<1e-8 hit or variant absent.",
        fontsize=7,
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    plt.close(fig)


def main() -> int:
    df = build_table()
    plot_heatmap(df)
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
