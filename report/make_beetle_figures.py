"""Figures and LaTeX tables for beetle_report.tex, from the beetle outputs. Run from claude_code_kit/
after beetle_check.py, beetle_model_compare.py, beetle_unscored.py and report/make_figures.py:

  python report/make_beetle_figures.py

Every count-area comparison uses the same 22 scored count areas with enough clear images.
"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, ".")
sys.path.insert(0, "report")
from beetle_rf import rank_r, score  # noqa: E402
from make_figures import FIGS, REGION_COLOUR, REGION_NAME, TABLES, style, write  # noqa: E402

B = Path("outputs/beetle_check")


def load():
    survey = pd.read_csv(B / "beetle_vs_survey.csv")
    check = pd.read_csv(B / "beetle_check.csv")
    preds = pd.read_csv(B / "model_comparison_predictions.csv")
    a = survey[survey.enough_images].reset_index(drop=True)
    a = a.merge(check[["Moor", "Count", "share_flagged__v0_original"] +
                      [c for c in check.columns if c.endswith("__v1_edges") and c.startswith("d")]],
                on=["Moor", "Count"])
    a = a.merge(preds[["Moor", "Count", "Gradient boosting (classic)", "main_combined"]], on=["Moor", "Count"])
    a["region_avg"] = [a[(a.GWCT_Region == r.GWCT_Region) & (a.Moor != r.Moor)].damage.mean()
                       for r in a.itertuples()]
    return survey, a, pd.read_csv(B / "unscored_estimates.csv")


STEPS = [("Browning only", "share_browned_raw"),
         ("Browning or greying, beyond\nthe moor-wide change", "share_flagged__v0_original"),
         ("+ 20 m margin round burns and\ncuts (rule-based estimate)", "satellite_estimate"),
         ("+ averaged with gradient\nboosting (final estimate)", "main_combined")]


def fig_progress(a):
    vals = [score(a, a[c].values)["typical_error"] for _, c in STEPS]
    bench = score(a, a.region_avg.values)["typical_error"]
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    y = np.arange(len(STEPS))[::-1]
    ax.barh(y, vals, height=0.58, color=["#9aa5b1", "#9aa5b1", "#5b9fd0", "#c1272d"], zorder=3)
    for yi, v in zip(y, vals):
        ax.text(v + 0.004, yi, f"{v:.2f}", va="center", fontsize=10)
    ax.axvline(bench, color="#444", ls="--", lw=1.1, zorder=4)
    ax.text(bench + 0.004, y[-1], f"region average,\nno satellite data ({bench:.2f})", fontsize=8.5,
            color="#444", va="center")
    ax.set_yticks(y, [s for s, _ in STEPS], fontsize=9.5)
    ax.set_xlabel("Typical difference from the survey (lower is better)")
    ax.set_xlim(0, max(vals) * 1.18)
    ax.grid(axis="x", alpha=0.25, zorder=0)
    fig.tight_layout()
    fig.savefig(FIGS / "beetle_fig_progress.png", dpi=220)
    plt.close(fig)


def fig_by_area(a, un):
    order = {"NYM": 0, "S.DALES/PEAK": 1, "N.DALES": 2}
    d = a.assign(o=a.GWCT_Region.map(order)).sort_values(["o", "damage", "main_combined"], ascending=[True, False, False])
    groups = [("NYM", "North York Moors"), ("S.DALES/PEAK", "Southern Dales and Peak"),
              ("N.DALES", "Northern Dales")]
    rows = []                                        # (label, survey, estimate, region); header rows have region None
    for code, label in groups:
        rows.append((label, np.nan, np.nan, None))
        rows += [(f"{r.Moor}, {r.Count}", r.damage, r.main_combined, code) for r in d[d.GWCT_Region == code].itertuples()]
    rows.append(("Not surveyed (Eggleston)", np.nan, np.nan, None))
    rows += [(f"{r.Moor}, {r.Count}", np.nan, r.combined, "unscored")
             for r in un.sort_values("combined", ascending=False).itertuples()]
    fig, ax = plt.subplots(figsize=(7.6, 0.26 * len(rows) + 0.9))
    y = np.arange(len(rows))[::-1]
    for yi, (name, s_, e, reg) in zip(y, rows):
        if reg is None:
            if yi != y[0]:
                ax.axhline(yi + 0.5, color="#bbb", lw=0.8)
            continue
        col = REGION_COLOUR.get(reg, "#7b8794")
        if np.isfinite(s_):
            ax.plot([s_, e], [yi, yi], color="#c9d1da", lw=2.2, zorder=1)
            ax.scatter(s_, yi, s=46, facecolor="white", edgecolor="#222", lw=1.3, zorder=3)
            ax.scatter(e, yi, s=46, color=col, edgecolor="white", lw=0.8, zorder=4)
        else:
            ax.scatter(e, yi, s=52, marker="D", color=col, edgecolor="white", lw=0.8, zorder=4)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=8.5)
    for t, r in zip(ax.get_yticklabels(), rows):
        if r[3] is None:
            t.set_fontweight("bold")
            t.set_fontstyle("italic")
            t.set_color("#333")
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(y[-1] - 0.7, y[0] + 0.6)
    ax.set_xlim(-0.03, 1.0)
    ax.set_xlabel("Proportion of count area damaged, summer 2025")
    ax.grid(axis="x", alpha=0.25, zorder=0)
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="#222", ms=7,
                      label="Survey, July 2025"),
               Line2D([], [], marker="o", ls="", color="#2b6cb0", ms=7, label="Satellite estimate"),
               Line2D([], [], marker="D", ls="", color="#7b8794", ms=6.5, label="Satellite estimate, not surveyed")]
    ax.legend(handles=handles, frameon=True, framealpha=1, edgecolor="#ddd", fontsize=8.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGS / "beetle_fig_by_area.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def tables(survey, a, un):
    n = json.loads(Path("outputs/report/report_numbers.json").read_text())

    # Survey summary by region
    reg = {"N.DALES": "Northern Dales", "NYM": "North York Moors", "S.DALES/PEAK": "Southern Dales and Peak"}
    rows = []
    for code in ["NYM", "S.DALES/PEAK", "N.DALES"]:
        g = survey[survey.GWCT_Region == code]
        rows.append([reg[code], len(g), f"{g.damage.mean():.2f}", f"{g.damage.min():.1f} to {g.damage.max():.1f}",
                     int((g.damage == 0).sum())])
    rows += ["mid", [r"\textbf{All}", len(survey), f"{survey.damage.mean():.2f}",
                     f"{survey.damage.min():.1f} to {survey.damage.max():.1f}", int((survey.damage == 0).sum())]]
    write("beetle_tab_survey", rows, "lrrrr",
          r"Region & \makecell[r]{Count areas\\scored} & \makecell[r]{Mean\\damage} & Range & \makecell[r]{Scored 0\\(no damage)}")

    # From first attempt to final estimate
    rows = []
    for label, col in STEPS + [("Gradient boosting alone", "Gradient boosting (classic)"),
                                ("Region average, no satellite data", "region_avg")]:
        s = score(a, a[col].values)
        lab = label.replace("\n", " ")
        if col == "main_combined":
            rows.append([rf"\textbf{{{lab}}}", rf"\textbf{{{s['typical_error']:.2f}}}", rf"\textbf{{{s['within_0.1']}}}",
                         rf"\textbf{{{s['within_0.2']}}}", rf"\textbf{{{s['rank_r']:.2f}}}"])
            rows.append("mid")
        else:
            rows.append([lab, f"{s['typical_error']:.2f}", s["within_0.1"], s["within_0.2"], f"{s['rank_r']:.2f}"])
    write("beetle_tab_progress", rows, ">{\\raggedright\\arraybackslash}p{7.2cm}rrrr",
          r"Estimate & \makecell[r]{Typical\\difference} & \makecell[r]{Within\\0.1} & \makecell[r]{Within\\0.2} & "
          r"\makecell[r]{Rank\\agreement}")

    # Which changes follow the survey (count-area medians, rank correlation)
    meas = [("Browning: NBR, summer 2025 against summer 2024", "dNBR_summer25_vs_24__v1_edges"),
            ("Greying: colour saturation, late against early summer 2025", "dSat_late_vs_early_2025__v1_edges"),
            ("Colour saturation, summer 2025 against summer 2024", "dSat_summer25_vs_24__v1_edges"),
            ("Greenness: NDVI, summer 2025 against summer 2024", "dNDVI_summer25_vs_24__v1_edges"),
            ("A year later: NBR, summer 2026 against summer 2024", "dNBR_summer26_vs_24__v1_edges"),
            (r"\textbf{Share of count area flagged (rule-based estimate)}", "satellite_estimate")]
    rows = [[lab, f"{rank_r(a.damage, a[c]):+.2f}".replace("-", "$-$")] for lab, c in meas]
    write("beetle_tab_signals", rows, "lr", r"Measure (median over the count area) & \makecell[r]{Rank correlation\\with the survey}")

    # Per count area, including the five not surveyed
    order = {"NYM": 0, "S.DALES/PEAK": 1, "N.DALES": 2}
    d = a.assign(o=a.GWCT_Region.map(order)).sort_values(["o", "damage", "Moor"])
    rows = []
    for code, g in d.groupby("o", sort=True):
        if rows:
            rows.append("mid")
        rows.append([rf"\multicolumn{{6}}{{l}}{{\textit{{{reg[g.GWCT_Region.iloc[0]]}}}}}"])
        for r in g.itertuples():
            gb = a.loc[r.Index, "Gradient boosting (classic)"]
            rows.append([f"{r.Moor}, {r.Count}", f"{r.damage:.1f}", f"{r.satellite_estimate:.2f}", f"{gb:.2f}",
                         rf"\textbf{{{r.main_combined:.2f}}}",
                         f"{r.main_combined - r.damage:+.2f}".replace("-", "$-$")])
    rows += ["mid", [r"\multicolumn{6}{l}{\textit{Not surveyed (Northern Dales)}}"]]
    for r in un.sort_values("combined", ascending=False).itertuples():
        rows.append([f"{r.Moor}, {r.Count}", "", f"{r.rule:.2f}", f"{r.gradient_boosting:.2f}",
                     rf"\textbf{{{r.combined:.2f}}}", ""])
    write("beetle_tab_areas", rows, "lrrrrr",
          r"Count area & Survey & Rule & \makecell[r]{Gradient\\boosting} & \makecell[r]{Final\\estimate} & Difference")

    # Beetle damage and management (from report_numbers.json)
    bm = n["beetle_management"]
    rows = [[s, f"{v['r_level']:+.2f}".replace("-", "$-$"), f"{v['r_vs_usual']:+.2f}".replace("-", "$-$")]
            for s, v in bm.items()]
    write("beetle_tab_management", rows, "lrr",
          r"Season & \makecell[r]{Rank correlation with\\\% burnt or cut} & \makecell[r]{Rank correlation with\\change from usual level}")
    print(f"Wrote beetle tables to {TABLES}")


if __name__ == "__main__":
    style()
    survey, a, un = load()
    fig_progress(a)
    fig_by_area(a, un)
    tables(survey, a, un)
    print(f"Wrote beetle figures to {FIGS}")
