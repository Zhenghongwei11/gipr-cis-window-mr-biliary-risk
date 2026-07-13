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
DOMAIN_SUMMARY = Path("docs/source_data/Source_Data_15_domain_summary.tsv")
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
    del df
    summary = pd.read_csv(DOMAIN_SUMMARY, sep="\t")
    summary["n_hits_p_lt_1e5"] = pd.to_numeric(summary["n_hits_p_lt_1e5"], errors="coerce")
    summary = summary[summary["rsid"].isin(VARIANTS)].copy()
    summary["variant_label"] = summary["rsid"].map({
        "rs10407429": "rs10407429\n(strict retained)",
        "rs17561351": "rs17561351\n(strict excluded)",
    })

    domain_order = (
        summary.groupby("domain")["n_hits_p_lt_1e5"]
        .sum()
        .sort_values()
        .index
        .tolist()
    )
    variants = ["rs10407429", "rs17561351"]
    colours = {"rs10407429": "#0072B2", "rs17561351": "#D55E00"}

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    y_base = list(range(len(domain_order)))
    offsets = {"rs10407429": 0.18, "rs17561351": -0.18}
    height = 0.32

    for variant in variants:
        sub = summary[summary["rsid"] == variant].set_index("domain")
        values = [float(sub.loc[d, "n_hits_p_lt_1e5"]) if d in sub.index else 0.0 for d in domain_order]
        y = [i + offsets[variant] for i in y_base]
        ax.barh(y, values, height=height, color=colours[variant], alpha=0.88, label=summary[summary["rsid"] == variant]["variant_label"].iloc[0])
        for yi, val in zip(y, values):
            if val > 0:
                ax.text(val + max(summary["n_hits_p_lt_1e5"]) * 0.015, yi, f"{int(val)}", va="center", fontsize=8)

    ax.set_yticks(y_base)
    ax.set_yticklabels(domain_order, fontsize=9)
    ax.set_xlabel("Retrieved association records (p < 1e-5)", fontsize=10)
    ax.set_ylabel("")
    ax.grid(axis="x", color="#DDDDDD", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01,
        0.01,
        "Counts are database association records after domain classification, not independent phenotype counts.",
        fontsize=7.5,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 1))

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    plt.close(fig)


def main() -> int:
    required_wald = [WALD_DIR / filename for _domain, _label, filename in WALD_TRAITS]
    if all(path.exists() for path in required_wald):
        df = build_table()
    elif OUT_TABLE.exists():
        df = pd.read_csv(OUT_TABLE, sep="\t")
    else:
        missing = ", ".join(str(path) for path in required_wald if not path.exists())
        raise FileNotFoundError(
            "Missing Wald diagnostic inputs and no cached association map is available: "
            f"{missing}"
        )
    plot_heatmap(df)
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
