#!/usr/bin/env python3
"""
Random Forest for heather beetle damage, labelled automatically from the July 2025 survey.

No hand-drawn labels: the survey supplies them. Every count area scored 0 is healthy ground;
count areas scored 0.5 or more are mostly damaged. The forest learns, pixel by pixel, what
separates them, from many Sentinel-2 measures instead of the two hand-picked rules in
beetle_check.py.

Step 1, Earth Engine (cached in rf_samples.csv; --resample to redo): up to --pixels random
pixels per count area, moorland only, excluding ground burnt or cut in 2024/25 or 2025/26
plus a 20 m margin. For each pixel, the change in all ten 10 m and 20 m bands (including
red-edge) and in NBR, NDVI, NDMI, colour saturation and brightness, for three comparisons
(summer 2025 vs 2024, late vs early summer 2025, summer 2026 vs 2024), plus the 2024 levels.

Step 2, local (scikit-learn), always leave-one-moor-out: when a moor is predicted, none of its
pixels or survey scores are used in training.
  regression   each pixel's target is its count area's survey score; a count area's
               estimate is the mean pixel prediction
  classifier   pixels from survey-0 count areas = healthy, survey >= 0.5 = damaged (others
               left out); a count area's estimate is its mean probability, mapped to a
               proportion with a straight line fitted on the training count areas
  combined     the mean of the classifier's estimate and the unfitted share flagged from
               beetle_check.py. They fail differently (the forest learns that bog greying is
               not damage; the rule keeps the full range on the greyed sites), so the mean
               does better than either. Chosen after seeing both, on 22 count areas.
Feature sets: "change" only (main), and "change + 2024 levels". Levels describe vegetation
type, so they let the forest learn region instead of damage; comparing the two shows whether
it does. The region shortcut is also checked by ranking within N.DALES alone.

Usage
  python beetle_rf.py                 # uses cached samples if present
  python beetle_rf.py --resample      # re-sample pixels in Earth Engine first
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

from beetle_check import COMPARISONS, WINDOWS
from burn_pipeline import DEFAULT_PROJECT, STRIP, Engine, ee_init, load_polygons, log, to_fc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BANDS10 = ["B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"]
INDICES = ["NBR", "NDVI", "NDMI", "Sat", "Bright"]
BUFFER_M = 20
CHUNK = 5                                           # count areas per Earth Engine request
REGION_COLOUR = {"N.DALES": "#2a78d6", "S.DALES/PEAK": "#eb6834", "NYM": "#1baf7a"}


# ---------------------------------------------------------------- step 1: sample pixels
def composite(ee, eng, start, end):
    med = eng.s2(start, end, bands=BANDS10).median()
    vis = med.select(["B2", "B3", "B4"])
    vmax = vis.reduce(ee.Reducer.max())
    return ee.Image.cat([med.select(BANDS10 + ["NBR", "NDVI"]),
                         med.normalizedDifference(["B8", "B11"]).rename("NDMI"),
                         vmax.subtract(vis.reduce(ee.Reducer.min())).divide(vmax).rename("Sat"),
                         vis.reduce(ee.Reducer.mean()).rename("Bright")])


def sample_pixels(ee, eng, uids, n_pixels):
    comp = {k: composite(ee, eng, *w) for k, w in WINDOWS.items()}
    names = BANDS10 + INDICES
    feats = [comp[a].subtract(comp[b]).rename([f"{n}__{c}" for n in names])
             for c, (a, b) in COMPARISONS.items()]
    feats.append(comp["late_2024"].rename([f"{n}__level2024" for n in names]))
    managed = eng.season(2024)["managed"].Or(eng.season(2025)["managed"])
    keep = eng.moorland.And(managed.focalMax(radius=BUFFER_M, units="meters").Not())
    img = ee.Image.cat(feats).updateMask(keep)

    def per_area(f):
        return (img.sample(region=f.geometry(), scale=10, numPixels=n_pixels, seed=7,
                           geometries=False, tileScale=4)
                .map(lambda p: p.set("uid", f.get("uid"))))

    rows = []
    for i in range(0, len(uids), CHUNK):
        part = eng.fc.filter(ee.Filter.inList("uid", uids[i:i + CHUNK]))
        got = part.map(per_area).flatten().getInfo()["features"]
        rows += [g["properties"] for g in got]
        log(f"  sampled count areas {i + 1}-{min(i + CHUNK, len(uids))} of {len(uids)}: {len(got)} pixels")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- step 2: models
def rank_r(a, b):
    return float(np.corrcoef(pd.Series(a).rank(), pd.Series(b).rank())[0, 1])


def area_means(px, values):
    return pd.Series(values, index=px.index).groupby([px.Moor, px.Count]).mean()


def lomo(px, areas, cols, kind):
    """Leave-one-moor-out estimates for every count area in `areas`."""
    est = {}
    for moor in areas.Moor.unique():
        tr, te = px[px.Moor != moor], px[px.Moor == moor]
        if kind == "regression":
            m = RandomForestRegressor(n_estimators=300, min_samples_leaf=20, max_features="sqrt",
                                      n_jobs=-1, random_state=0).fit(tr[cols], tr.damage)
            est.update(area_means(te, m.predict(te[cols])).to_dict())
        else:
            lab = tr[(tr.damage == 0) | (tr.damage >= 0.5)]
            m = RandomForestClassifier(n_estimators=300, min_samples_leaf=20, max_features="sqrt",
                                       oob_score=True, n_jobs=-1, random_state=0)
            m.fit(lab[cols], lab.damage >= 0.5)
            # training count areas: out-of-bag probability where labelled, direct prediction otherwise
            p_tr = pd.Series(m.predict_proba(tr[cols])[:, 1], index=tr.index)
            p_tr.loc[lab.index] = m.oob_decision_function_[:, 1]
            cal = area_means(tr, p_tr.values).to_frame("p").join(
                tr.groupby(["Moor", "Count"]).damage.first())
            slope, icpt = np.polyfit(cal.p, cal.damage, 1)
            p_te = area_means(te, m.predict_proba(te[cols])[:, 1])
            est.update((icpt + slope * p_te).to_dict())
    return areas.set_index(["Moor", "Count"]).index.map(est).to_numpy(dtype=float).clip(0, 1)


def score(areas, est):
    e = np.abs(est - areas.damage.values)
    nd = (areas.GWCT_Region == "N.DALES").values
    return {"typical_error": round(e.mean(), 3), "within_0.2": int((e <= 0.2001).sum()),
            "within_0.1": int((e <= 0.1001).sum()), "rank_r": round(rank_r(areas.damage, est), 2),
            "rank_r_N.DALES_only": round(rank_r(areas.damage[nd], est[nd]), 2)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--polygons", default="data/count_areas_v6.gpkg")
    ap.add_argument("--project", default=DEFAULT_PROJECT)
    ap.add_argument("--threshold", type=float, default=0.17, help="for the managed mask")
    ap.add_argument("--pixels", type=int, default=300, help="pixels sampled per count area")
    ap.add_argument("--resample", action="store_true")
    ap.add_argument("--out", default="outputs/beetle_check")
    a = ap.parse_args()
    out = Path(a.out)
    cache = out / "rf_samples.csv"

    gdf_all = load_polygons(a.polygons)
    core = gdf_all[gdf_all.Count_Type != STRIP]
    if a.resample or not cache.exists():
        ee = ee_init(a.project)
        eng = Engine(ee, to_fc(ee, core), a.threshold, 0.1, 40, False, aoi_fc=to_fc(ee, gdf_all))
        log(f"Sampling up to {a.pixels} pixels in each of {len(core)} count areas ...")
        sample_pixels(ee, eng, [int(u) for u in core.uid], a.pixels).to_csv(cache, index=False)
    px = pd.read_csv(cache).merge(core[["uid", "Moor", "Count", "GWCT_Region"]], on="uid")

    survey = pd.read_csv(out / "beetle_vs_survey.csv")
    areas = survey[survey.enough_images].reset_index(drop=True)
    px = px.merge(areas[["Moor", "Count", "damage"]], on=["Moor", "Count"])
    px = px.dropna().reset_index(drop=True)
    log(f"{len(px):,} pixels from {px.groupby(['Moor', 'Count']).ngroups} scored count areas "
        f"(median {px.groupby(['Moor', 'Count']).size().median():.0f} per area)")

    change = [c for c in px.columns if any(c.endswith(f"__{k}") for k in COMPARISONS)]
    level = [c for c in px.columns if c.endswith("__level2024")]
    runs = {"RF regression, change only": (change, "regression"),
            "RF classifier, change only": (change, "classifier"),
            "RF regression, change + 2024 levels": (change + level, "regression"),
            "RF classifier, change + 2024 levels": (change + level, "classifier")}
    rows = [{"method": "Current method (unfitted share flagged)", **score(areas, areas.satellite_estimate.values)}]
    preds = areas[["GWCT_Region", "Moor", "Count", "damage", "satellite_estimate"]].copy()
    for name, (cols, kind) in runs.items():
        est = lomo(px, areas, cols, kind)
        preds[name] = est.round(2)
        rows.append({"method": name, **score(areas, est)})
        log(f"  {name}: done")
    combined = (areas.satellite_estimate.values + preds["RF classifier, change only"].values) / 2
    preds["Combined (current + RF classifier)"] = combined.round(2)
    rows.append({"method": "Combined (current + RF classifier)", **score(areas, combined)})
    res = pd.DataFrame(rows)
    res.to_csv(out / "rf_results.csv", index=False)
    preds.to_csv(out / "rf_predictions.csv", index=False)
    log("\n" + res.to_string(index=False))

    # Which measures matter: a regression forest on all scored count areas, change features only
    m = RandomForestRegressor(n_estimators=500, min_samples_leaf=20, max_features="sqrt",
                              n_jobs=-1, random_state=0).fit(px[change], px.damage)
    imp = pd.Series(m.feature_importances_, index=change).sort_values(ascending=False)
    imp.round(4).to_csv(out / "rf_importance.csv", header=["importance"])
    log("\nTop measures:\n" + imp.head(12).round(3).to_string())

    # Figure: current method against the combined estimate, plus importances
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    res_i = res.set_index("method")
    for ax, (col, method, title) in zip(axes[:2], [
            ("satellite_estimate", "Current method (unfitted share flagged)", "Rule-based estimate (nothing fitted)"),
            ("Combined (current + RF classifier)", "Combined (current + RF classifier)",
             "Rule-based estimate averaged with the Random Forest")]):
        r = res_i.loc[method]
        ax.plot([0, 1], [0, 1], color="#999", lw=0.9, zorder=1)
        ax.fill_between([0, 1], [-0.2, 0.8], [0.2, 1.2], color="#eee", zorder=0)
        for reg, g in preds.groupby("GWCT_Region"):
            ax.scatter(g.damage, g[col], s=40, color=REGION_COLOUR.get(reg, "#888"), edgecolor="white",
                       lw=1.0, label=reg, zorder=3)
        for _, row in preds.iterrows():
            ax.annotate(row.Count, (row.damage, row[col]), xytext=(4, 3), textcoords="offset points",
                        fontsize=6.5, color="#444")
        ax.set_title(f"{title}\ntypical error {r.typical_error:.2f}, within 0.2: {int(r['within_0.2'])}/{len(preds)}, "
                     f"within 0.1: {int(r['within_0.1'])}/{len(preds)}, rank r {r.rank_r:.2f}",
                     loc="left", fontsize=9.5, weight="bold")
        ax.set_xlabel("Surveyed beetle damage, July 2025")
        ax.set_xlim(-0.05, 1.0)
        ax.set_ylim(-0.05, 1.0)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Satellite estimate (each moor predicted without its own data)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    top = imp.head(12)[::-1]
    axes[2].barh([c.replace("__", ", ") for c in top.index], top.values, color="#2a78d6")
    axes[2].set_title("Measures the forest relies on most", loc="left", fontsize=9.5, weight="bold")
    axes[2].tick_params(labelsize=8)
    axes[2].spines[["top", "right"]].set_visible(False)
    fig.suptitle("Beetle damage: Random Forest labelled automatically from the survey, tested one moor at a time. "
                 "Grey band = within 0.2",
                 x=0.01, ha="left", fontsize=12, weight="bold")
    plt.tight_layout()
    fig.savefig(out / "fig_beetle_rf.png", dpi=160, facecolor="white")
    log(f"Wrote {out}")


if __name__ == "__main__":
    main()
