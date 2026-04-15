import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-chemy")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
THERMO = DATA / "thermo"
PRESENTATION = ROOT / "presentation"
FIGURES = PRESENTATION / "figures"
TABLES = PRESENTATION / "tables"

EXPERIMENTS = {
    "Elements only": THERMO / "experiments" / "elements_only",
    "All BURCAT": THERMO / "experiments" / "burcat_all",
    "80/20 split": THERMO / "experiments" / "burcat_split_80_20_seed42",
}


def ensure_dirs():
    FIGURES.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)


def load_json(path):
    with open(path) as f:
        return json.load(f)


def load_jsonl(path):
    if not path.exists():
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def latex_escape(value):
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def fmt(value, digits=1):
    if value is None or not finite(value):
        return "--"
    return f"{value:.{digits}f}"


def savefig(name):
    path = FIGURES / name
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    return path


def set_style():
    plt.rcParams.update(
        {
            "figure.figsize": (8.2, 4.6),
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "legend.fontsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "savefig.bbox": "tight",
            "savefig.dpi": 180,
            "pdf.fonttype": 42,
        }
    )


def load_compound_names():
    names = {}
    for subdir in [DATA / "chems" / "mapped", DATA / "chems" / "unmapped"]:
        for path in sorted(subdir.glob("*.jsonl")):
            for entry in load_jsonl(path):
                name = entry.get("cmpdname") or entry.get("iupacname") or f"CID {entry['cid']}"
                names[entry["cid"]] = name
    return names


def make_pipeline_figure():
    labels = [
        "Compound\ncatalog",
        "LLM reaction\nprompts",
        "Validation\nand parsing",
        "Balanced\nreaction graph",
        "LLM reaction\nthermo",
        "BURCAT\nreferences",
        "Weighted\nformation solve",
        "Held-out\nvalidation",
    ]
    x = [0, 1.45, 2.9, 4.35, 5.8, 5.8, 7.25, 8.7]
    y = [0, 0, 0, 0, 0.7, -0.7, 0, 0]
    colors = ["#e9f2f7", "#e9f2f7", "#e9f2f7", "#e9f2f7", "#f7eddc", "#eee7f7", "#e6f4e8", "#f7e6e6"]

    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.axis("off")
    for i, (label, xi, yi, color) in enumerate(zip(labels, x, y, colors)):
        ax.text(
            xi,
            yi,
            label,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.06", facecolor=color, edgecolor="#3a3a3a"),
        )
        if i in {0, 1, 2, 3, 6}:
            ax.annotate("", xy=(x[i + 1] - 0.48, y[i + 1]), xytext=(xi + 0.48, yi), arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.annotate("", xy=(7.25 - 0.52, 0.05), xytext=(5.8 + 0.5, 0.6), arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.annotate("", xy=(7.25 - 0.52, -0.05), xytext=(5.8 + 0.5, -0.6), arrowprops=dict(arrowstyle="->", lw=1.4))
    ax.set_xlim(-0.7, 9.4)
    ax.set_ylim(-1.45, 1.35)
    savefig("pipeline.pdf")


def make_llm_evidence_figures():
    estimates = load_jsonl(THERMO / "evidence" / "reaction_estimates.jsonl")
    dH = np.array([e["dH"] for e in estimates if finite(e.get("dH")) and abs(e["dH"]) < 500])
    dG = np.array([e["dG"] for e in estimates if finite(e.get("dG")) and abs(e["dG"]) < 500])

    plt.figure(figsize=(8.2, 4.4))
    bins = np.linspace(-250, 250, 80)
    plt.hist(dH, bins=bins, alpha=0.58, label=r"$\Delta H_r$", color="#2f6f9f")
    plt.hist(dG, bins=bins, alpha=0.58, label=r"$\Delta G_r$", color="#c65b3a")
    plt.xlabel("LLM reaction estimate, kcal per balanced equation")
    plt.ylabel("Reaction count")
    plt.title("Distribution of LLM reaction thermochemistry")
    plt.legend()
    savefig("llm_reaction_distribution.pdf")

    means_h = np.array([e["dH"] for e in estimates if finite(e.get("dH")) and finite(e.get("std_dH")) and abs(e["dH"]) < 500])
    std_h = np.array([e["std_dH"] for e in estimates if finite(e.get("dH")) and finite(e.get("std_dH")) and abs(e["dH"]) < 500])
    means_g = np.array([e["dG"] for e in estimates if finite(e.get("dG")) and finite(e.get("std_dG")) and abs(e["dG"]) < 500])
    std_g = np.array([e["std_dG"] for e in estimates if finite(e.get("dG")) and finite(e.get("std_dG")) and abs(e["dG"]) < 500])

    plt.figure(figsize=(8.2, 4.4))
    plt.scatter(np.abs(means_h), std_h, s=8, alpha=0.22, label=r"$\Delta H_r$", color="#2f6f9f")
    plt.scatter(np.abs(means_g), std_g, s=8, alpha=0.22, label=r"$\Delta G_r$", color="#c65b3a")
    plt.xlabel("Absolute mean estimate, kcal")
    plt.ylabel("Across-sample standard deviation, kcal")
    plt.title("LLM uncertainty grows with reaction magnitude")
    plt.legend()
    savefig("llm_uncertainty.pdf")


def make_burcat_summary_figure():
    report = load_json(THERMO / "reports" / "burcat_parse_report.json")
    labels = ["Species\nparsed", "Parse\nfailures", "Element\nrefs", "CAS\nrefs", "Refs\nwritten"]
    values = [
        report["species_entries"],
        report["parse_failure_count"],
        report["element_references"],
        report["cas_references"],
        report["references_written"],
    ]
    plt.figure(figsize=(7.2, 4.2))
    bars = plt.bar(labels, values, color=["#557a95", "#c94c4c", "#88a85b", "#7d6ca7", "#5b9c8a"])
    plt.ylabel("Count")
    plt.title("BURCAT extraction coverage")
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, value, f"{value:,}", ha="center", va="bottom", fontsize=10)
    savefig("burcat_summary.pdf")


def metric(exp_dir, subset, prop, key):
    return load_json(exp_dir / "metrics.json")[subset][prop][key]


def make_metrics_figure():
    labels = ["Elements only\n(test)", "All BURCAT\n(overall)", "80/20 split\n(test)"]
    dirs = [EXPERIMENTS["Elements only"], EXPERIMENTS["All BURCAT"], EXPERIMENTS["80/20 split"]]
    subsets = ["test", "overall", "test"]
    x = np.arange(len(labels))
    width = 0.34
    h_mae = [metric(d, s, "dHf", "mae") for d, s in zip(dirs, subsets)]
    g_mae = [metric(d, s, "dGf", "mae") for d, s in zip(dirs, subsets)]

    plt.figure(figsize=(8.2, 4.4))
    plt.bar(x - width / 2, h_mae, width, label=r"$\Delta H_f$", color="#2f6f9f")
    plt.bar(x + width / 2, g_mae, width, label=r"$\Delta G_f$", color="#c65b3a")
    plt.xticks(x, labels)
    plt.ylabel("MAE vs BURCAT, kcal/mol")
    plt.title("Reference policy changes validation error")
    plt.legend()
    savefig("policy_mae.pdf")


def rows_for(exp, subset=None):
    rows = load_jsonl(exp / "validation.jsonl")
    if subset is not None:
        rows = [row for row in rows if row.get("subset") == subset]
    return rows


def make_scatter_and_error_figures():
    split = EXPERIMENTS["80/20 split"]
    rows = rows_for(split, "test")
    for prop, title, filename in [
        ("dHf", r"$\Delta H_f^\circ$", "dhf_scatter_split_test.pdf"),
        ("dGf", r"$\Delta G_f^\circ$", "dgf_scatter_split_test.pdf"),
    ]:
        xs = np.array([r[f"{prop}_ref"] for r in rows if finite(r.get(f"{prop}_ref")) and finite(r.get(f"{prop}_solved"))])
        ys = np.array([r[f"{prop}_solved"] for r in rows if finite(r.get(f"{prop}_ref")) and finite(r.get(f"{prop}_solved"))])
        plt.figure(figsize=(5.6, 5.1))
        plt.scatter(xs, ys, s=18, alpha=0.68, color="#2f6f9f", edgecolors="none")
        lo = float(min(xs.min(), ys.min()))
        hi = float(max(xs.max(), ys.max()))
        pad = 0.05 * (hi - lo)
        plt.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#444444", lw=1.2, linestyle="--")
        plt.xlabel(f"BURCAT reference {title}, kcal/mol")
        plt.ylabel(f"Solved {title}, kcal/mol")
        plt.title(f"Held-out split: {title} prediction")
        plt.xlim(lo - pad, hi + pad)
        plt.ylim(lo - pad, hi + pad)
        savefig(filename)

    configs = [
        ("Elements only", rows_for(EXPERIMENTS["Elements only"], "test")),
        ("80/20 train", rows_for(split, "train")),
        ("80/20 test", rows_for(split, "test")),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4), sharey=False)
    for ax, prop, title in zip(axes, ["dHf", "dGf"], [r"$\Delta H_f$", r"$\Delta G_f$"]):
        data = []
        labels = []
        for label, source_rows in configs:
            errors = [r[f"{prop}_error"] for r in source_rows if finite(r.get(f"{prop}_error"))]
            data.append(errors)
            labels.append(label)
        ax.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True)
        ax.set_title(f"{title} error")
        ax.set_ylabel("Solved - reference, kcal/mol")
        ax.tick_params(axis="x", rotation=15)
        ax.grid(True, axis="y", alpha=0.25)
    savefig("error_distributions.pdf")


def make_reaction_residuals():
    estimates = {e["rid"]: e for e in load_jsonl(THERMO / "evidence" / "reaction_estimates.jsonl")}
    corrected = load_jsonl(EXPERIMENTS["80/20 split"] / "reactions_corrected.jsonl")
    h = []
    g = []
    for row in corrected:
        est = estimates.get(row["rid"])
        if not est:
            continue
        if finite(row.get("dH")) and finite(est.get("dH")):
            h.append(row["dH"] - est["dH"])
        if finite(row.get("dG")) and finite(est.get("dG")):
            g.append(row["dG"] - est["dG"])
    plt.figure(figsize=(8.2, 4.4))
    bins = np.linspace(-120, 120, 90)
    plt.hist(h, bins=bins, alpha=0.58, label=r"$\Delta H_r$", color="#2f6f9f")
    plt.hist(g, bins=bins, alpha=0.58, label=r"$\Delta G_r$", color="#c65b3a")
    plt.xlabel("Solved reaction value minus LLM estimate, kcal")
    plt.ylabel("Reaction count")
    plt.title("Reaction-level correction induced by the global solve")
    plt.legend()
    savefig("reaction_residuals.pdf")


def write_counts_table():
    refs = load_jsonl(THERMO / "evidence" / "formation_references.jsonl")
    estimates = load_jsonl(THERMO / "evidence" / "reaction_estimates.jsonl")
    report = load_json(THERMO / "reports" / "burcat_parse_report.json")
    split_config = load_json(EXPERIMENTS["80/20 split"] / "config.json")
    rows = [
        ("LLM reaction estimates", f"{len(estimates):,}"),
        ("Reactions used in solves", f"{split_config['reactions']:,}"),
        ("Compounds solved in 80/20 split", f"{load_json(EXPERIMENTS['80/20 split'] / 'metadata.json')['chemicals']:,}"),
        ("BURCAT species parsed", f"{report['species_entries']:,}"),
        ("BURCAT parse failures", f"{report['parse_failure_count']:,}"),
        ("Element references", f"{sum(1 for r in refs if r.get('reference_type') == 'element'):,}"),
        ("Non-element BURCAT references", f"{sum(1 for r in refs if r.get('reference_type') != 'element'):,}"),
        ("80/20 split train references", f"{split_config['formation_references_train']:,}"),
        ("80/20 split test references", f"{split_config['formation_references_test']:,}"),
    ]
    write_tabular(TABLES / "counts.tex", ["Quantity", "Value"], rows)


def write_metrics_table():
    rows = []
    specs = [
        ("Elements only", "test"),
        ("All BURCAT", "overall"),
        ("80/20 split train", "train"),
        ("80/20 split test", "test"),
    ]
    dir_by_label = {
        "Elements only": EXPERIMENTS["Elements only"],
        "All BURCAT": EXPERIMENTS["All BURCAT"],
        "80/20 split train": EXPERIMENTS["80/20 split"],
        "80/20 split test": EXPERIMENTS["80/20 split"],
    }
    for label, subset in specs:
        metrics = load_json(dir_by_label[label] / "metrics.json")[subset]
        rows.append(
            (
                label,
                f"{metrics['dHf']['count']}/{metrics['dHf']['reference_count']}",
                fmt(metrics["dHf"]["mae"]),
                fmt(metrics["dHf"]["rmse"]),
                f"{metrics['dGf']['count']}/{metrics['dGf']['reference_count']}",
                fmt(metrics["dGf"]["mae"]),
                fmt(metrics["dGf"]["rmse"]),
            )
        )
    write_tabular(
        TABLES / "metrics.tex",
        [
            "Policy / subset",
            r"$H$ n",
            r"$H$ MAE",
            r"$H$ RMSE",
            r"$G$ n",
            r"$G$ MAE",
            r"$G$ RMSE",
        ],
        rows,
    )


def write_examples_table():
    names = load_compound_names()
    rows = []
    for row in rows_for(EXPERIMENTS["80/20 split"], "test"):
        if not finite(row.get("dHf_error")) or not finite(row.get("dGf_error")):
            continue
        score = abs(row["dHf_error"]) + abs(row["dGf_error"])
        rows.append((score, row))
    rows = sorted(rows, key=lambda item: item[0])[:8]
    table_rows = []
    for _score, row in rows:
        table_rows.append(
            (
                names.get(row["cid"], f"CID {row['cid']}")[:32],
                row.get("burcat_species") or f"CID {row['cid']}",
                fmt(row["dHf_ref"]),
                fmt(row["dHf_solved"]),
                fmt(row["dHf_error"]),
                fmt(row["dGf_error"]),
            )
        )
    write_tabular(
        TABLES / "examples.tex",
        [
            "Compound",
            "BURCAT",
            r"$H_f$ ref",
            r"$H_f$ solved",
            r"$H_f$ err",
            r"$G_f$ err",
        ],
        table_rows,
    )


def write_tabular(path, headers, rows):
    spec = "l" + "r" * (len(headers) - 1)
    lines = [rf"\begin{{tabular}}{{{spec}}}", r"\toprule"]
    lines.append(" & ".join(latex_escape(h) if not str(h).startswith("$") else str(h) for h in headers) + r" \\")
    lines.append(r"\midrule")
    for row in rows:
        lines.append(" & ".join(latex_escape(cell) for cell in row) + r" \\")
    lines.extend([r"\bottomrule", r"\end{tabular}", ""])
    path.write_text("\n".join(lines))


def main():
    ensure_dirs()
    set_style()
    make_pipeline_figure()
    make_llm_evidence_figures()
    make_burcat_summary_figure()
    make_metrics_figure()
    make_scatter_and_error_figures()
    make_reaction_residuals()
    write_counts_table()
    write_metrics_table()
    write_examples_table()
    print(f"Wrote figures to {FIGURES}")
    print(f"Wrote tables to {TABLES}")


if __name__ == "__main__":
    main()
