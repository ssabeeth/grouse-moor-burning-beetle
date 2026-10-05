#!/usr/bin/env python3
"""
Every number the v6 report quotes, computed from the run outputs, plus the threshold
sensitivity figure. Run after the pipeline runs (main, t013, t015, t019, no_s2c) and
beetle_check.py. Writes outputs/report/report_numbers.json, a few CSVs and a figure.

Usage
  python report_numbers.py
"""

import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from burn_pipeline import REGION, STRIP

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RUNS = {0.13: "outputs_t013", 0.15: "outputs_t015", 0.17: "outputs", 0.19: "outputs_t019"}
MAIN = 0.17
OUT = Path("outputs/report")
UNCHANGED = [("Dallowgill", "Bishops"), ("Ramsgill", "Raygill"), ("Swinton", "Potts"), ("Grinton", "Long Gill"),
             ("Snilesworth", "Lodge"), ("Danby", "Low Moor"), ("Moor House", "Behind House")]
SHEWRING = {"North York Moors": 6.3, "North Pennines": 2.7}          # %/yr, burning only (Shewring et al. 2024)
THRESH_COLOUR = {0.13: "#a6cbe3", 0.15: "#5b9fd0", 0.17: "#2171b5", 0.19: "#08306b"}   # light to dark = low to high


def load(folder):
    df = pd.read_csv(Path(folder) / "results_by_season.csv")
    df = df[df.Count_Type != STRIP].copy()
    df["Region"] = df.Moor.map(REGION).fillna("Other")
    return df


def pct(g):
    return 100 * g.managed_ha.sum() / g.area_ha.sum()


def rank_r(a, b):
    return round(float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1]), 2)


def region_only_mae(d):
    """Benchmark without satellite data: predict each moor's damage as the mean of the other
    moors in its GWCT region (leave-one-moor-out)."""
    err = []
    for r in d.itertuples():
        same = d[(d.GWCT_Region == r.GWCT_Region) & (d.Moor != r.Moor)]
        rest = same if len(same) else d[d.Moor != r.Moor]
        err.append(abs(rest.damage.mean() - r.damage))
    return round(float(np.mean(err)), 2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    R = {t: load(f) for t, f in RUNS.items()}
    m = R[MAIN]
    areas = m.drop_duplicates(["Moor", "Count"])
    n = {"count_areas": len(areas), "moors": areas.Moor.nunique(), "area_ha": round(areas.area_ha.sum()),
         "seasons": m.season.nunique(), "first_season": m.season.min(), "last_season": m.season.max()}

    # Seasons
    S = m.groupby("season").agg(managed_ha=("managed_ha", "sum"), area_ha=("area_ha", "sum"),
                                images=("post_images_per_pixel", "median"))
    S["pct"] = 100 * S.managed_ha / S.area_ha
    n["mean_ha_per_yr"] = round(S.managed_ha.mean())
    n["mean_pct_per_yr"] = round(pct(m), 1)
    n["season_max"], n["season_min"] = S.pct.idxmax(), S.pct.idxmin()
    pre = S[S.index < "2021"]
    post = S[(S.index >= "2021") & (S.index < "2024")]
    n["pre2021_mean_ha"], n["post2021_mean_ha"] = round(pre.managed_ha.mean()), round(post.managed_ha.mean())

    # Sentinel-2C
    s2c = load("outputs_no_s2c")
    with_c, without_c = S.loc["2024/25", "managed_ha"], s2c.managed_ha.sum()
    n["s2c"] = {"with_ha": round(with_c), "without_ha": round(without_c),
                "change_pct": round(100 * (without_c / with_c - 1)),
                "next_highest_ha": round(S.drop("2024/25").managed_ha.max()),
                "rank_r": rank_r(m[m.pre_year == 2024].set_index(["Moor", "Count"]).managed_ha,
                                 s2c.set_index(["Moor", "Count"]).managed_ha.reindex(
                                     m[m.pre_year == 2024].set_index(["Moor", "Count"]).index))}
    S.round(2).to_csv(OUT / "seasons.csv")

    # Moors, regions, controls
    M = m.groupby(["Moor", "Region"]).apply(pct).rename("pct").reset_index().sort_values("pct", ascending=False)
    M.round(2).to_csv(OUT / "moors.csv", index=False)
    G = m.groupby("Region").apply(pct)
    n["regions"] = {r: round(v, 1) for r, v in G.items()}
    n["shewring"] = SHEWRING
    n["nym_np_ratio"] = round(G["North York Moors"] / G["North Pennines"], 2)
    n["controls"] = {k: round(M.set_index("Moor").pct[k], 2) for k in
                     ["Moor House", "Geltsdale", "Crossgill", "Wemmergill"]}

    # Thresholds
    T = pd.DataFrame({t: df.groupby("season").managed_ha.sum() for t, df in R.items()})
    TP = 100 * T / S.area_ha.values[:, None]
    TM = pd.DataFrame({t: df.groupby("Moor").apply(pct) for t, df in R.items()})
    TG = pd.DataFrame({t: df.groupby("Region").apply(pct) for t, df in R.items()})
    n["thresholds"] = {str(t): {
        "total_change_pct": round(100 * (T[t].sum() / T[MAIN].sum() - 1)),
        "season_rank_r": rank_r(T[t], T[MAIN]), "moor_rank_r": rank_r(TM[t], TM[MAIN]),
        "nym_np_ratio": round(TG.loc["North York Moors", t] / TG.loc["North Pennines", t], 2),
        "post_pre_ratio": round(T[t][T.index.isin(post.index)].mean() / T[t][T.index.isin(pre.index)].mean(), 2),
        "moor_house": round(TM.loc["Moor House", t], 2), "geltsdale": round(TM.loc["Geltsdale", t], 2)}
        for t in RUNS}
    TP.round(2).to_csv(OUT / "threshold_seasons_pct.csv")

    # Eggleston boundary
    b = pd.read_csv("outputs/summary_boundary_sensitivity.csv").set_index(["scope", "boundary"]).mean_annual_pct
    n["eggleston"] = {"count_areas": round(b[("Eggleston", "count areas")], 1),
                      "whole_block": round(b[("Eggleston", "count areas + strip")], 1),
                      "np_count_areas": round(b[("North Pennines", "count areas")], 1),
                      "np_whole_block": round(b[("North Pennines", "count areas + strip")], 1)}

    # v5
    v5 = pd.read_csv("reference/results_v5.csv")
    j = m[m.pre_year <= 2024].merge(v5[["Moor", "Count", "season", "area_ha", "burn_ha", "cut_wf_ha", "total_ha"]],
                                    on=["Moor", "Count", "season"], suffixes=("", "_v5"))
    u = j[[(a, c) in UNCHANGED for a, c in zip(j.Moor, j.Count)]]
    n["v5"] = {"shared_count_areas": j[["Moor", "Count"]].drop_duplicates().shape[0],
               "pct_v5": round(100 * j.total_ha.sum() / j.area_ha_v5.sum(), 2),
               "pct_v6": round(100 * j.managed_ha.sum() / j.area_ha.sum(), 2),
               "unchanged_red_ha": round(u.burn_ha.sum()), "unchanged_blue_ha": round(u.cut_wf_ha.sum()),
               "unchanged_v6_ha": round(u.managed_ha.sum()),
               "zero_share_v5": round(100 * (j.total_ha < 0.05).mean()),
               "zero_share_v6": round(100 * (j.managed_ha < 0.05).mean())}
    w = m[m.Count == "Waskerley"].set_index("season").pct_of_area.round(1)
    n["waskerley"] = w.to_dict()
    n["pikestone_2022"] = round(m[(m.Count == "Pikestone Fell") & (m.pre_year == 2022)].pct_of_area.iloc[0], 1)
    n["grinton_2024_images"] = float(m[(m.Count == "Long Gill") & (m.pre_year == 2024)].post_images_per_pixel.iloc[0])
    n["median_images"] = float(m.post_images_per_pixel.median())

    # Beetle: management vs survey, and satellite vs survey
    bs = pd.read_csv("outputs/beetle_check/beetle_vs_survey.csv")
    x = m.merge(bs[["Moor", "Count", "damage"]], on=["Moor", "Count"])
    usual = x[x.pre_year < 2024].groupby(["Moor", "Count"]).pct_of_area.mean()
    bl = {}
    for yr in (2023, 2024, 2025):
        s = x[x.pre_year == yr].set_index(["Moor", "Count"])
        bl[s.season.iloc[0]] = {"r_level": rank_r(s.damage, s.pct_of_area),
                                "r_vs_usual": rank_r(s.damage, s.pct_of_area - usual.reindex(s.index))}
    n["beetle_management"] = bl
    models = pd.read_csv("outputs/beetle_check/beetle_models.csv")
    ok = bs[bs.enough_images]
    err = (ok.satellite_estimate - ok.damage).abs()
    n["beetle_survey"] = {
        "n_scored": len(bs), "n_good_images": len(ok), "mae": round(err.mean(), 2),
        "within_0.2": int((err <= 0.2001).sum()), "within_0.1": int((err <= 0.1001).sum()),
        "rank_r": rank_r(ok.damage, ok.satellite_estimate),
        "raw_mae": round((ok.share_browned_raw - ok.damage).abs().mean(), 2),   # same 22 areas as the rest
        "raw_rank_r": rank_r(ok.damage, ok.share_browned_raw),
        "low_damage_bias": round((ok.satellite_estimate - ok.damage)[ok.damage <= 0.1].mean(), 2),
        "region_only_mae": region_only_mae(ok),
        "models": models.to_dict("records")}

    rf = pd.read_csv("outputs/beetle_check/rf_results.csv").set_index("method")
    imp = pd.read_csv("outputs/beetle_check/rf_importance.csv", index_col=0).importance
    n["beetle_rf"] = {k: rf.loc[m].to_dict() for k, m in [
        ("rule", "Current method (unfitted share flagged)"), ("rf_classifier", "RF classifier, change only"),
        ("rf_regression", "RF regression, change only"), ("rf_with_levels", "RF classifier, change + 2024 levels"),
        ("combined", "Combined (current + RF classifier)")]}
    n["beetle_rf"]["red_edge_share_of_top10"] = int(sum(c.split("__")[0] in ("B5", "B6", "B7") for c in imp.index[:10]))
    rp = pd.read_csv("outputs/beetle_check/rf_predictions.csv").set_index("Count")
    n["beetle_rf"]["trout_beck"] = {c: float(rp.loc["Trout Beck", c]) for c in
                                    ["satellite_estimate", "RF classifier, change only", "Combined (current + RF classifier)"]}

    gp = pd.read_csv("outputs/beetle_check/model_comparison_predictions.csv").set_index("Count")
    gi = pd.read_csv("outputs/beetle_check/gb_importance.csv", index_col=0).importance
    n["beetle_gb"] = {"trout_beck_alone": float(gp.loc["Trout Beck", "Gradient boosting (classic)"]),
                      "trout_beck_combined": float(gp.loc["Trout Beck", "main_combined"]),
                      "low_damage_bias": round(float((gp.main_combined - gp.damage)[gp.damage <= 0.1].mean()), 3),
                      "red_edge_in_top10": int(sum(c.split("__")[0] in ("B5", "B6", "B7") for c in gi.index[:10])),
                      "region_holdout": pd.read_csv("outputs/beetle_check/region_holdout.csv").to_dict("records"),
                      "nested": pd.read_csv("outputs/beetle_check/nested_selection.csv").to_dict("records"),
                      "nested_choices": pd.read_csv("outputs/beetle_check/nested_selection_choices.csv")
                      .model_chosen.value_counts().to_dict()}
    mc = pd.read_csv("outputs/beetle_check/model_comparison.csv")
    n["beetle_models"] = mc[["model", "combined_with_rule", "typical_error", "within_0.1", "vs_rule_90pc"]].to_dict("records")

    (OUT / "report_numbers.json").write_text(json.dumps(n, indent=2, default=str))
    print(json.dumps({k: v for k, v in n.items() if k not in ("beetle_survey",)}, indent=1, default=str))
    print(json.dumps({k: v for k, v in n["beetle_survey"].items() if k != "models"}, default=str))

    # Figure: season totals at each threshold
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(9, 4.6))
    x_pos = np.arange(len(TP))
    for t in sorted(RUNS, reverse=True):
        main = t == MAIN
        ax.plot(x_pos, TP[t], marker="o", ms=5 if main else 4, lw=2.4 if main else 1.6,
                color=THRESH_COLOUR[t], zorder=3 if main else 2)
        ax.text(x_pos[-1] + 0.15, TP[t].iloc[-1], f"{t:.2f}" + (" (main)" if main else ""),
                va="center", fontsize=8.5, color="#333")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(TP.index, fontsize=9)
    ax.set_xlim(-0.3, len(TP) - 0.2)
    ax.set_ylim(0, TP.values.max() * 1.1)
    ax.set_ylabel("% of count area burnt or cut")
    ax.grid(axis="y", alpha=0.25)
    ax.set_title("Season totals at four thresholds: levels shift, the pattern does not",
                 loc="left", fontsize=11.5, weight="bold")
    plt.tight_layout()
    fig.savefig(OUT / "fig_threshold_sensitivity.png", dpi=200)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
