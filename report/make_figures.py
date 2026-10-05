"""Caption-ready figures and LaTeX tables for report_v6.tex. Figures have no in-figure titles
(the captions carry them); tables are written from the run outputs, so no number is typed by hand.

Run from claude_code_kit/ after the pipeline runs, report_numbers.py, beetle_model_compare.py
and v5_to_v6_changes.py --summarise:

  python report/make_figures.py

Writes report/figures/ and report/tables/. Image sheets (contact sheets, evidence sheets) are copied as JPEG.
"""
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, ".")
import burn_pipeline as bp  # noqa: E402

FIGS = Path("report/figures")
INK, MUTED = "#222222", "#666666"
REGION_COLOUR = {"N.DALES": "#2b6cb0", "NYM": "#2a9d72", "S.DALES/PEAK": "#e0663a"}
MENTIONED = ["Low Moor", "Ever Rigg", "Raygill", "Bishops", "Trout Beck", "South Side", "Waskerley"]
REGION_NAME = {"N.DALES": "Northern Dales", "NYM": "North York Moors",
               "S.DALES/PEAK": "Southern Dales and Peak"}


def place_labels(ax, xs, ys, names, fs=7):
    """Label each point without overlaps: try spots right, left, above and below, then further
    out with a thin leader line. Measured with the renderer, so it works at any figure size."""
    fig = ax.figure
    fig.tight_layout()                                          # place labels on the final layout
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    pts = ax.transData.transform(np.column_stack([xs, ys]))
    r = 3.8 * fig.dpi / 72                                      # marker radius in pixels
    placed = []
    spots = [(6, 0, "left"), (-6, 0, "right")] + [(dx, dy, ha) for k in (1, 2, 3, 4)
             for dy in (11 * k, -11 * k) for dx, ha in ((6, "left"), (-6, "right"))]
    inside = ax.get_window_extent(rend)
    for i in np.argsort(-np.asarray(ys)):
        if not names[i]:
            continue
        for dx, dy, ha in spots:
            t = ax.annotate(str(names[i]), (xs[i], ys[i]), xytext=(dx, dy), textcoords="offset points",
                            ha=ha, va="center", fontsize=fs)
            bb = t.get_window_extent(rend).expanded(1.06, 1.15)
            t.remove()
            out = bb.x0 < inside.x0 or bb.x1 > inside.x1 or bb.y0 < inside.y0 or bb.y1 > inside.y1
            hit = out or any(bb.overlaps(b) for b in placed) or any(
                np.hypot(min(max(px, bb.x0), bb.x1) - px, min(max(py, bb.y0), bb.y1) - py) < r
                for px, py in pts)
            if not hit:
                break
        ax.annotate(str(names[i]), (xs[i], ys[i]), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va="center", fontsize=fs, color="#444",
                    arrowprops=dict(arrowstyle="-", color="#aaa", lw=0.5, shrinkA=0, shrinkB=3) if dy else None)
        placed.append(bb)


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": "#444",
                         "axes.labelcolor": INK, "xtick.color": "#444", "ytick.color": "#444",
                         "axes.titlelocation": "left", "axes.titleweight": "bold",
                         "axes.titlesize": 10.5})


def main_charts():
    df = pd.read_csv("outputs/results_by_season.csv")
    gdf = bp.load_polygons("data/count_areas_v6.gpkg", None, None)
    bp.make_charts(df, gdf, Path("outputs"), figs=FIGS, titles=False)


def thresholds():
    tp = pd.read_csv("outputs/report/threshold_seasons_pct.csv", index_col="season")
    colour = {"0.13": "#a6cbe3", "0.15": "#5b9fd0", "0.17": "#2171b5", "0.19": "#08306b"}
    fig, ax = plt.subplots(figsize=(9, 4.2))
    x = np.arange(len(tp))
    for t in sorted(tp.columns, reverse=True):
        main = t == "0.17"
        ax.plot(x, tp[t], marker="o", ms=5 if main else 4, lw=2.6 if main else 1.6,
                color=colour[t], zorder=3 if main else 2)
        ax.text(x[-1] + 0.15, tp[t].iloc[-1], f"{t}" + (" (main)" if main else ""),
                va="center", fontsize=9, color=INK)
    ax.set_xticks(x, tp.index, fontsize=9.5)
    ax.set_xlim(-0.3, len(tp) - 0.1)
    ax.set_ylim(0, tp.values.max() * 1.1)
    ax.set_ylabel("% of count area burnt or cut")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_threshold_sensitivity.png", dpi=220)
    plt.close(fig)


def beetle():
    p = pd.read_csv("outputs/beetle_check/model_comparison_predictions.csv")
    imp = pd.read_csv("outputs/beetle_check/gb_importance.csv", index_col=0).importance
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.9))
    for ax, (col, title) in zip(axes, [("rule", "(a) Rule-based estimate"),
                                       ("main_combined", "(b) Averaged with gradient boosting")]):
        ax.plot([0, 1], [0, 1], color="#999", lw=0.9, zorder=1)
        ax.fill_between([0, 1], [-0.1, 0.9], [0.1, 1.1], color="#e9eef5", zorder=0)
        for region, g in p.groupby("GWCT_Region"):
            ax.scatter(g.damage, g[col], s=42, color=REGION_COLOUR[region], edgecolor="white",
                       lw=1.0, label=REGION_NAME[region], zorder=3)
        ax.set_title(title)
        ax.set_xlabel("Surveyed beetle damage, July 2025")
        ax.set_xlim(-0.05, 1.0)
        ax.set_ylim(-0.05, 1.0)
        ax.set_aspect("equal")
    axes[0].set_ylabel("Satellite estimate")
    axes[0].legend(frameon=False, fontsize=8.5, loc="upper left")
    for ax, col in zip(axes, ["rule", "main_combined"]):
        named = (p[col] - p.damage).abs().gt(0.1) | p.Count.isin(MENTIONED)   # selective labels
        place_labels(ax, p.damage.values, p[col].values, p.Count.where(named, "").values, fs=8)
    fig.savefig(FIGS / "fig_beetle_gb.png", dpi=220, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    # Which measures gradient boosting relies on
    top = imp.head(10)[::-1]
    red_edge = [c.split("__")[0] in ("B5", "B6", "B7") for c in top.index]
    nice = {"late_vs_early_2025": "late vs early summer 2025", "summer25_vs_24": "summer 2025 vs 2024",
            "summer26_vs_24": "summer 2026 vs 2024"}
    fig, ax = plt.subplots(figsize=(7.2, 3.7))
    ax.barh([f"{c.split('__')[0]}, {nice[c.split('__')[1]]}" for c in top.index],
            top.values, color=["#c1272d" if r else "#7b8794" for r in red_edge])
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#c1272d", label="red-edge band (B5 to B7)"),
                       Patch(color="#7b8794", label="other band or index")],
              frameon=False, fontsize=9, loc="lower right")
    ax.set_xlabel("Importance (share of total)")
    ax.tick_params(axis="y", labelsize=9.5)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_beetle_importance.png", dpi=220, facecolor="white")
    plt.close(fig)


def v5_to_v6():
    t = pd.read_csv("outputs/v5_to_v6/v5_to_v6_summary.csv").set_index("run")
    s = pd.read_csv("outputs/v5_to_v6/v5_to_v6_by_season.csv", index_col="season")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4.0), gridspec_kw={"width_ratios": [0.8, 1.45]})
    steps = [("v1 method,\nv1 bound-\naries", "v5 boundaries, v5 settings", "#9aa5b1"),
             ("v1 method,\nv2 bound-\naries", "v6 boundaries, v5 settings", "#9aa5b1"),
             ("+ window\nfrom\n16 April", "v6 boundaries, 16 April window only", "#5b9fd0"),
             ("+ 0.1 ha\npatch\n(v2)", "v6 boundaries, v6 settings (main run)", "#c1272d")]
    vals = [t.loc[k, "pct_per_year"] for _, k, _ in steps]
    xs = np.arange(len(steps))
    a1.bar(xs, vals, width=0.62, color=[c for *_, c in steps], zorder=3)
    for x, v in zip(xs, vals):
        a1.text(x, v + max(vals) * 0.02, f"{v:.1f}%", ha="center", va="bottom", fontsize=9.5)
    a1.set_xticks(xs, [lbl for lbl, *_ in steps], fontsize=8.5)
    a1.set_ylabel("% of count area burnt or cut per season")
    a1.set_ylim(0, max(vals) * 1.18)
    a1.grid(axis="y", alpha=0.25, zorder=0)
    a1.set_title("(a) One change at a time")

    x = np.arange(len(s))
    w = 0.27
    cols = [("1 Apr window", "#9aa5b1", "v1 method: from 1 April, 0.5 ha"),
            ("16 Apr window", "#5b9fd0", "from 16 April, 0.5 ha"),
            ("16 Apr + 0.1 ha (v6)", "#c1272d", "v2: from 16 April, 0.1 ha")]
    for i, (c, colour, label) in enumerate(cols):
        a2.bar(x + (i - 1) * w, s[c], width=w, color=colour, label=label, zorder=3)
    a2.set_xticks(x, s.index, fontsize=8.5, rotation=45, ha="right")
    a2.set_ylabel("Burnt or cut (ha)")
    a2.grid(axis="y", alpha=0.25, zorder=0)
    a2.legend(frameon=False, fontsize=8.5, loc="upper left")
    a2.set_ylim(0, s.values.max() * 1.3)
    a2.set_title("(b) By season, v2 boundaries")
    fig.tight_layout()
    fig.savefig(FIGS / "fig_v5_to_v6.png", dpi=220)
    plt.close(fig)

def image_sheets():
    """Satellite image sheets as high-quality JPEG: a quarter of the PNG size, so the PDF can be emailed."""
    from PIL import Image
    src = list(Path("outputs_images/figures/contact_sheets").glob("*.png")) + \
        list(Path("outputs/beetle_check").glob("fig_evidence_*.png"))
    for f in src:
        Image.open(f).convert("RGB").save(FIGS / f"{f.stem}.jpg", quality=88, optimize=True)


# ---------------------------------------------------------------- tables
TABLES = Path("report/tables")


def fmt_ha(v):
    return f"{v:,.0f}"


def write(name, rows, cols, head, note=None):
    """A booktabs tabular body: `cols` column spec, `head` header row, rows as lists of cells."""
    lines = [rf"\begin{{tabular}}{{{cols}}}", r"\toprule", head + r" \\", r"\midrule"]
    for r in rows:
        if r == "mid":
            lines.append(r"\midrule")
        else:
            lines.append(" & ".join(str(c) for c in r) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (TABLES / f"{name}.tex").write_text("\n".join(lines) + "\n")


def tables():
    import json
    TABLES.mkdir(parents=True, exist_ok=True)
    n = json.loads(Path("outputs/report/report_numbers.json").read_text())

    # Seasons
    se = pd.read_csv("outputs/report/seasons.csv")
    rows = [[r.season, fmt_ha(r.managed_ha), f"{100 * r.managed_ha / r.area_ha:.1f}", f"{r.images:.1f}"]
            for r in se.itertuples()]
    rows += ["mid", [r"\textbf{Mean}", rf"\textbf{{{se.managed_ha.mean():.0f}}}",
                     rf"\textbf{{{100 * se.managed_ha.sum() / se.area_ha.sum():.1f}}}", f"{n['median_images']:.1f}"]]
    write("tab_seasons", rows, "lrrr",
          r"Season & Burnt or cut (ha) & \% of count area & Clear images per pixel\textsuperscript{a}")

    # Moors
    d = pd.read_csv("outputs/results_by_season.csv")
    d = d[d.Count_Type != "Strip"]
    m = d.groupby("Moor").agg(n=("Count", "nunique"), h=("managed_ha", "sum"), a_seasons=("area_ha", "sum"))
    m["area"] = d.groupby(["Moor", "Count"]).area_ha.first().groupby("Moor").sum()
    m["pct"] = 100 * m.h / m.a_seasons
    m["Region"] = m.index.map(bp.REGION)
    m = m.sort_values("pct", ascending=False)
    rows = [[moor, r.Region.replace("Dales / Nidderdale", "Dales and Nidderdale"), int(r.n), fmt_ha(r.area),
             f"{r.pct:.1f}" if r.pct >= 1 else f"{r.pct:.2f}"] for moor, r in m.iterrows()]
    write("tab_moors", rows, "llrrr", r"Moor & Region & Count areas & Area (ha) & Mean \% a year")

    # Thresholds (appendix)
    th = n["thresholds"]
    rows = []
    for t in ["0.13", "0.15", "0.17", "0.19"]:
        v = th[t]
        main = t == "0.17"
        rows.append([rf"\textbf{{{t}}}" if main else t,
                     "main run" if main else f"{v['total_change_pct']:+d}\\%".replace("-", "$-$"),
                     f"{v['season_rank_r']:.2f}", f"{v['moor_rank_r']:.2f}", f"{v['nym_np_ratio']:.2f}",
                     f"{v['post_pre_ratio']:.2f}", f"{v['moor_house']:.2f}", f"{v['geltsdale']:.1f}"])
    write("tab_thresholds", rows, "lrrrrrrr",
          r"Threshold & \makecell[r]{Total\\vs 0.17} & \makecell[r]{Season\\order\textsuperscript{a}} & "
          r"\makecell[r]{Moor\\order\textsuperscript{a}} & \makecell[r]{NYM /\\NP} & \makecell[r]{After /\\before 2021} & "
          r"\makecell[r]{Moor\\House \%} & \makecell[r]{Geltsdale\\\%}")

    # v5 -> v6 steps
    t = pd.read_csv("outputs/v5_to_v6/v5_to_v6_summary.csv")
    label = {"v5 as reported (burn class only)": "Version 1 as published (main class)",
             "v5 boundaries, v5 settings": "Version 1 method, rebuilt here",
             "v5 boundaries, v6 settings": "Version 2 method",
             "v6 boundaries, v5 settings": "Version 1 method",
             "v6 boundaries, 16 April start only": r"\quad start on 16 April only",
             "v6 boundaries, 16 April window only": r"\quad start on 16 April, end on 30 June",
             "v6 boundaries, 0.1 ha patch only": r"\quad 0.1 ha patch only",
             "v6 boundaries, v6 settings (main run)": r"\textbf{Version 2}"}
    t = t.set_index("run")
    groups = [(r"\textit{Version 1 boundaries (2,546 ha)}", ["v5 as reported (burn class only)",
               "v5 boundaries, v5 settings", "v5 boundaries, v6 settings"]),
              (r"\textit{Version 2 boundaries (4,218 ha)}", ["v6 boundaries, v5 settings",
               "v6 boundaries, 16 April start only", "v6 boundaries, 16 April window only",
               "v6 boundaries, 0.1 ha patch only",
               "v6 boundaries, v6 settings (main run)"])]
    rows = []
    for head, keys in groups:
        if rows:
            rows.append("mid")
        rows.append([rf"\multicolumn{{7}}{{l}}{{{head}}}"])
        for k in keys:
            r = t.loc[k]
            win = {"04-01,06-01": "1 Apr to 31 May", "04-16,06-01": "16 Apr to 31 May",
                   "04-16,07-01": "16 Apr to 30 Jun"}[r.after_window]
            rows.append([label[k], win, f"{r.min_patch_ha:.1f}", f"{r.pct_per_year:.1f}",
                         f"{r.pct_seasons_zero:.0f}", fmt_ha(r.unchanged7_total_ha), fmt_ha(r.s2024_ha)])
    write("tab_v5_to_v6", rows, "llrrrrr",
          r"Run & After window & \makecell[r]{Min.\\patch (ha)} & \makecell[r]{\%\\a year} & "
          r"\makecell[r]{Seasons with\\none (\%)} & \makecell[r]{Seven\\areas (ha)\textsuperscript{a}} & "
          r"\makecell[r]{2024/25\\(ha)}")

    # Beetle models (appendix)
    mc = pd.read_csv("outputs/beetle_check/model_comparison.csv")
    rows = []
    rule = mc[mc.model.str.startswith("Rule")].iloc[0]
    rows.append(["Rule-based estimate only", f"{rule.typical_error:.3f}", int(rule["within_0.1"]),
                 f"{rule.rank_r:.2f}", "", "", "", "", ""])
    rows.append("mid")
    names = {"Gradient boosting (classic)": r"\textbf{Gradient boosting (classic)}",
             "Gradient boosting (histogram, LightGBM-style)": "Gradient boosting (LightGBM style)",
             "SVM (RBF kernel)": "Support vector machine", "Random Forest": "Random forest",
             "Logistic regression": "Logistic regression", "Extra Trees": "Extra trees"}
    comb = mc[mc.combined_with_rule].set_index("model")
    alone = mc[~mc.combined_with_rule & ~mc.model.str.startswith("Rule")].set_index("model")
    for k in comb.sort_values("typical_error").index:
        a, c = alone.loc[k], comb.loc[k]
        ci = c.vs_rule_90pc.replace(" to ", " to ").replace("-", "$-$").replace("+", "$+$")
        rows.append([names[k], f"{a.typical_error:.3f}", int(a["within_0.1"]), f"{a.rank_r:.2f}",
                     f"{c.typical_error:.3f}", int(c["within_0.1"]), int(c["within_0.2"]), f"{c.rank_r:.2f}", ci])
    write("tab_beetle_models", rows, "lrrrrrrrl",
          r"& \multicolumn{3}{c}{Model alone} & \multicolumn{5}{c}{Averaged with the rule} \\ \cmidrule(lr){2-4}\cmidrule(lr){5-9}"
          r"Model & \makecell[r]{Typical\\diff.} & \makecell[r]{Within\\0.1} & \makecell[r]{Rank\\$\rho$} & "
          r"\makecell[r]{Typical\\diff.} & \makecell[r]{Within\\0.1} & \makecell[r]{Within\\0.2} & \makecell[r]{Rank\\$\rho$} "
          r"& \makecell[l]{90\% interval\\vs rule}")

    # Beetle per count area (appendix)
    pr = pd.read_csv("outputs/beetle_check/model_comparison_predictions.csv")
    reg = {"N.DALES": "Northern Dales", "NYM": "North York Moors", "S.DALES/PEAK": "Southern Dales and Peak"}
    pr = pr.sort_values(["GWCT_Region", "damage", "Moor"])
    rows = []
    for region, g in pr.groupby("GWCT_Region", sort=False):
        if rows:
            rows.append("mid")
        rows.append([rf"\multicolumn{{6}}{{l}}{{\textit{{{reg[region]}}}}}"])
        for _, r in g.iterrows():
            gb = r["Gradient boosting (classic)"]
            rows.append([f"{r.Moor}, {r.Count}", f"{r.damage:.1f}", f"{r.rule:.2f}", f"{gb:.2f}",
                         rf"\textbf{{{r.main_combined:.2f}}}", f"{r.main_combined - r.damage:+.2f}".replace("-", "$-$")])
    write("tab_beetle_areas", rows, "lrrrrr",
          r"Count area & Survey & Rule & \makecell[r]{Gradient\\boosting} & Combined & Difference")

    # Eggleston boundary comparison (appendix)
    eb = pd.read_csv("outputs/summary_boundary_sensitivity.csv")
    rows = [[r.scope, r.boundary.replace("count areas + strip", "count areas + strip (Eleanor's block)"),
             fmt_ha(r.area_ha), f"{r.mean_annual_pct:.1f}"] for r in eb.itertuples()]
    write("tab_eggleston", rows, "llrr", r"Scope & Boundary & Area (ha) & Mean \% a year")
    print(f"Wrote {TABLES}")


if __name__ == "__main__":
    FIGS.mkdir(parents=True, exist_ok=True)
    style()
    main_charts()
    style()
    thresholds()
    beetle()
    image_sheets()
    if Path("outputs/v5_to_v6/v5_to_v6_summary.csv").exists():
        v5_to_v6()
    tables()
    print(f"Wrote {FIGS}")
