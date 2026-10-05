#!/usr/bin/env python3
"""
Can Sentinel-2 see heather beetle damage, and how close does it get to the survey numbers?

Compares each count area's change over summer 2025 and 2026 with the July 2025 beetle
damage survey (reference/HEATHER_BEETLE_DAMAGE_SUMMER_2025.xlsx). Three comparisons:

  summer25_vs_24       late summer (Jul to Sep) 2025 minus late summer 2024
  late_vs_early_2025   late summer 2025 minus early summer (16 Apr to 30 Jun) 2025
  summer26_vs_24       late summer 2026 minus late summer 2024

Beetle-killed heather goes red-brown in the summer it is attacked, then grey, so four measures:

  dNBR, dNDVI   negative = lost greenness and moisture (browning or greying)
  dBright       mean visible reflectance; positive = paler
  dSat          colour saturation of the visible bands; negative = flatter colour (greyer)

Each measure is computed under three masks, to show what each fix contributes:

  v0_original        moorland, minus ground burnt or cut in 2024/25 or 2025/26
  v1_edges           as v0, also minus a 20 m margin around that ground (burn and cut edges)
  v2_edges_evergreen as v1, heather-like ground only: the more evergreen half of moorland pixels
                     (smallest summer-minus-winter NDVI swing, 2022 and 2023, before the outbreak).
                     Heather stays dark all year; bog grasses and cotton-grass green up and bleach,
                     and their seasonal colour change is not beetle damage. (ESA WorldCover cannot
                     do this: it labels 98.5% of the count areas grassland and none shrubland.)

Fix 3: count areas with fewer than 3 clear images per pixel in any window are flagged and
left out of the "filtered" comparison.

Survey comparison, two estimates per variant:
  share flagged   (headline) share of usable pixels that browned (NBR, summer 2025 vs 2024) by
                  more than 0.10 beyond the moor-wide change, or greyed (saturation, late vs
                  early summer 2025) by more than 0.05 beyond it. Read directly as a proportion;
                  nothing is fitted to the survey.
  calibrated      straight-line model on the same two measures, each moor predicted from a fit
                  to the other moors only (leave-one-moor-out).
"before" in the figure is the first attempt: raw share of pixels whose NBR fell by more than 0.10,
with no moor-wide correction.

Usage
  python beetle_check.py --polygons data/count_areas_v6.gpkg
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from burn_pipeline import DEFAULT_PROJECT, STRIP, Engine, ee_init, load_polygons, log, to_fc

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BEETLE = Path(__file__).parent / "reference/HEATHER_BEETLE_DAMAGE_SUMMER_2025.xlsx"
NAME_FIX = {"greencastle": "greencastles"}          # beetle sheet spelling -> count area name
REGION_COLOUR = {"N.DALES": "#2a78d6", "S.DALES/PEAK": "#eb6834", "NYM": "#1baf7a"}

WINDOWS = {"late_2024": ("2024-07-01", "2024-10-01"),
           "early_2025": ("2025-04-16", "2025-07-01"),
           "late_2025": ("2025-07-01", "2025-10-01"),
           "late_2026": ("2026-07-01", "2026-10-01")}
COMPARISONS = {"summer25_vs_24": ("late_2025", "late_2024"),
               "late_vs_early_2025": ("late_2025", "early_2025"),
               "summer26_vs_24": ("late_2026", "late_2024")}
METRICS = {"dNBR": "NBR", "dNDVI": "NDVI", "dBright": "Bright", "dSat": "Sat"}

VARIANTS = {"v0_original": dict(buffer_m=0, evergreen=False),
            "v1_edges": dict(buffer_m=20, evergreen=False),
            "v2_edges_evergreen": dict(buffer_m=20, evergreen=True)}
SUMMER = [("2022-06-15", "2022-09-01"), ("2023-06-15", "2023-09-01")]   # pre-outbreak seasonal swing
WINTER = [("2022-02-15", "2022-04-10"), ("2023-02-15", "2023-04-10")]
MIN_CLEAR = 3                                       # clear images per pixel, per window (fix 3)
MODEL = ["dNBR_summer25_vs_24", "dSat_late_vs_early_2025"]
FLAG = {"dNBR_summer25_vs_24": -0.10, "dSat_late_vs_early_2025": -0.05}   # beyond the moor-wide change


def key(s):
    k = str(s).lower().replace(" ", "")
    return NAME_FIX.get(k, k)


def rank_corr(x, y, n_perm=20000, seed=0):
    """Spearman correlation and a two-sided permutation p-value (no scipy needed)."""
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    r = np.corrcoef(rx, ry)[0, 1]
    rng = np.random.default_rng(seed)
    null = np.array([np.corrcoef(rx, rng.permutation(ry))[0, 1] for _ in range(n_perm)])
    return r, (np.abs(null) >= abs(r) - 1e-12).mean()


def composite(ee, eng, start, end):
    """Median NBR, NDVI, visible brightness and visible colour saturation."""
    med = eng.s2(start, end).select(["B2", "B3", "B4", "NBR", "NDVI"]).median()
    vis = med.select(["B2", "B3", "B4"])
    vmax = vis.reduce(ee.Reducer.max())
    sat = vmax.subtract(vis.reduce(ee.Reducer.min())).divide(vmax)   # 0 = grey, higher = coloured
    return ee.Image.cat([med.select(["NBR", "NDVI"]), vis.reduce(ee.Reducer.mean()).rename("Bright"),
                         sat.rename("Sat")])


def evergreen_mask(ee, eng):
    """Heather-like pixels: summer-minus-winter NDVI swing below the median across the count areas."""
    def ndvi(windows):
        col = eng.s2(*windows[0])
        for w in windows[1:]:
            col = col.merge(eng.s2(*w))
        return col.select("NDVI").median()
    swing = ndvi(SUMMER).subtract(ndvi(WINTER)).updateMask(eng.moorland).rename("swing")
    cut = swing.reduceRegion(reducer=ee.Reducer.median(), geometry=eng.aoi, scale=20,
                             maxPixels=1e9, bestEffort=True, tileScale=4).get("swing")
    return swing.lt(ee.Number(cut)).unmask(0)


def measure(ee, eng):
    managed = eng.season(2024)["managed"].Or(eng.season(2025)["managed"])   # burnt or cut 2024/25, 2025/26
    evergreen = evergreen_mask(ee, eng)
    comp = {k: composite(ee, eng, *w) for k, w in WINDOWS.items()}
    names, diffs = [], []
    for name, (a, b) in COMPARISONS.items():
        for metric, band in METRICS.items():
            names.append(f"{metric}_{name}")
            diffs.append(comp[a].select(band).subtract(comp[b].select(band)))
    diff = ee.Image.cat(diffs).rename(names)

    bands = [eng.s2(*w).select("NBR").count().unmask(0).rename(f"clear_{k}") for k, w in WINDOWS.items()]
    bands.append(diff.select("dNBR_summer25_vs_24").lt(-0.10)
                 .updateMask(eng.moorland.And(managed.Not())).rename("share_browned_raw"))
    for v, cfg in VARIANTS.items():
        drop = managed.focalMax(radius=cfg["buffer_m"], units="meters") if cfg["buffer_m"] else managed
        keep = eng.moorland.And(evergreen if cfg["evergreen"] else 1).And(drop.Not())
        masked = diff.updateMask(keep)
        off = masked.select(list(FLAG)).reduceRegion(reducer=ee.Reducer.median(), geometry=eng.aoi,
                                                     scale=20, maxPixels=1e9, bestEffort=True, tileScale=4)
        flagged = None
        for col, thr in FLAG.items():
            f = masked.select(col).subtract(ee.Number(off.get(col))).lt(thr)
            flagged = f if flagged is None else flagged.Or(f)
        bands += [masked.rename([f"{n}__{v}" for n in names]),
                  flagged.rename(f"share_flagged__{v}"),
                  keep.unmask(0).rename(f"share_used__{v}")]
    img = ee.Image.cat(bands)
    fc = img.reduceRegions(collection=eng.fc, reducer=ee.Reducer.mean(), scale=10, tileScale=8).getInfo()
    return pd.DataFrame([f["properties"] for f in fc["features"]])


def lomo(d, features):
    """Leave-one-moor-out straight-line predictions of surveyed damage, clipped to 0..1."""
    y = d.damage.values
    X = np.column_stack([np.ones(len(d))] + [d[f].values for f in features])
    pred = np.zeros(len(d))
    for moor in d.Moor.unique():
        te = (d.Moor == moor).values
        beta, *_ = np.linalg.lstsq(X[~te], y[~te], rcond=None)
        pred[te] = X[te] @ beta
    return np.clip(pred, 0, 1)


def scores(y, p):
    e = np.abs(p - y)
    return {"n": len(y), "MAE": round(e.mean(), 3), "bias": round((p - y).mean(), 3),
            "within_0.1": round((e <= 0.1001).mean(), 2), "within_0.2": round((e <= 0.2001).mean(), 2),
            "spearman": round(rank_corr(y, p, n_perm=1)[0], 2)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--polygons", default="data/count_areas_v6.gpkg")
    ap.add_argument("--project", default=DEFAULT_PROJECT)
    ap.add_argument("--threshold", type=float, default=0.17, help="for the managed masks")
    ap.add_argument("--out", default="outputs/beetle_check")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    gdf = load_polygons(a.polygons)
    ee = ee_init(a.project)
    eng = Engine(ee, to_fc(ee, gdf), a.threshold, 0.1, 40, False)
    log("Measuring summer 2024, 2025 and 2026 change per count area, three mask variants (one request) ...")
    sat = measure(ee, eng)
    sat = sat.merge(gdf[["uid", "Count_Type", "GWCT_Region", "area_ha"]], on="uid")
    sat = sat[sat.Count_Type != STRIP]

    beetle = pd.read_excel(BEETLE)
    beetle = beetle.rename(columns={beetle.columns[-1]: "damage"})
    beetle["k"] = beetle.MOOR.map(key) + "|" + beetle.COUNT.map(key)
    sat["k"] = sat.Moor.map(key) + "|" + sat.Count.map(key)
    if unmatched := set(beetle.k) - set(sat.k):
        log(f"WARNING: beetle rows with no count area: {sorted(unmatched)}")
    df = sat.merge(beetle[["k", "damage"]], on="k", how="left").drop(columns="k")
    clear_cols = ["clear_late_2024", "clear_early_2025", "clear_late_2025"]       # windows the model uses
    df["min_clear"] = df[clear_cols].min(axis=1).round(1)
    df["enough_images"] = df.min_clear >= MIN_CLEAR
    df.sort_values(["GWCT_Region", "Moor", "Count"]).to_csv(out / "beetle_check.csv", index=False)

    scored = df.dropna(subset=["damage"]).reset_index(drop=True)
    log(f"{len(scored)} count areas with a beetle score. Too few clear images (< {MIN_CLEAR}): "
        f"{', '.join(scored[~scored.enough_images].Count) or 'none'}")
    for v in VARIANTS:
        used = scored[f"share_used__{v}"]
        log(f"  {v}: share of count area used, median {used.median():.0%}, lowest "
            f"{scored.loc[used.idxmin(), 'Count']} {used.min():.0%}")

    # Correlations, per variant
    rows = []
    for v in VARIANTS:
        for col in [*MODEL, "dNBR_summer26_vs_24", "dSat_summer25_vs_24", "share_flagged"]:
            c = f"{col}__{v}"
            sub = scored.dropna(subset=[c])
            r, p = rank_corr(sub.damage, sub[c])
            rows.append({"variant": v, "measure": col, "n": len(sub), "spearman_r": round(r, 2), "p_perm": round(p, 4)})
    pd.DataFrame(rows).to_csv(out / "beetle_correlations.csv", index=False)

    # Survey comparison: calibrated (leave-one-moor-out) and uncalibrated, all sites and filtered
    rows, preds = [], {}
    for v in VARIANTS:
        for subset, sub in [("all", scored), ("enough images", scored[scored.enough_images])]:
            sub = sub.dropna(subset=[f"{f}__{v}" for f in MODEL]).reset_index(drop=True)
            p = lomo(sub, [f"{f}__{v}" for f in MODEL])
            rows.append({"variant": v, "sites": subset, "estimate": "calibrated (leave-one-moor-out)",
                         **scores(sub.damage.values, p)})
            rows.append({"variant": v, "sites": subset, "estimate": "uncalibrated share flagged",
                         **scores(sub.damage.values, sub[f"share_flagged__{v}"].values)})
            if subset == "all":
                preds[v] = pd.Series(p, index=pd.MultiIndex.from_frame(sub[["Moor", "Count"]]))
    res = pd.DataFrame(rows)
    res.to_csv(out / "beetle_models.csv", index=False)
    log("\n" + res.to_string(index=False))

    final = "v1_edges"
    comp = scored[["GWCT_Region", "Moor", "Count", "damage", "min_clear", "enough_images",
                   f"share_used__{final}"]].copy()
    for v in VARIANTS:
        comp[f"estimate_{v}"] = comp.set_index(["Moor", "Count"]).index.map(preds[v]).values.round(2)
    comp["share_browned_raw"] = scored.share_browned_raw.round(2)
    comp["satellite_estimate"] = scored[f"share_flagged__{final}"].round(2)        # headline
    comp["error"] = (comp.satellite_estimate - comp.damage).round(2)
    comp.sort_values("damage").to_csv(out / "beetle_vs_survey.csv", index=False)

    # Figure: before (raw share browned) and after (share flagged, final variant)
    raw = scores(comp.damage.values, comp.share_browned_raw.values)
    ok = comp.enough_images
    fin = scores(comp.damage[ok].values, comp.satellite_estimate[ok].values)
    panels = [("share_browned_raw", f"Before: raw share browned\ntypical error {raw['MAE']:.2f}, "
                                    f"within 0.2: {raw['within_0.2']:.0%}, rank r {raw['spearman']:.2f}"),
              ("satellite_estimate", f"After: browned or greyed beyond moor-wide change, burn edges out\n"
                                     f"typical error {fin['MAE']:.2f}, within 0.2: {fin['within_0.2']:.0%}, "
                                     f"rank r {fin['spearman']:.2f} (enough-image sites)")]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.4), sharex=True, sharey=True)
    for ax, (col, title) in zip(axes, panels):
        ax.plot([0, 1], [0, 1], color="#999", lw=0.9, zorder=1)
        ax.fill_between([0, 1], [-0.2, 0.8], [0.2, 1.2], color="#eee", zorder=0)
        for reg, g in comp.groupby("GWCT_Region"):
            good = g.enough_images
            ax.scatter(g.damage[good], g[col][good], s=42, color=REGION_COLOUR.get(reg, "#888"),
                       edgecolor="white", lw=1.0, label=reg, zorder=3)
            ax.scatter(g.damage[~good], g[col][~good], s=42, facecolor="none",
                       edgecolor=REGION_COLOUR.get(reg, "#888"), lw=1.3, zorder=3)
        for r in comp.itertuples():
            ax.annotate(r.Count, (r.damage, getattr(r, col)), xytext=(4, 3),
                        textcoords="offset points", fontsize=6.5, color="#444")
        ax.set_title(title, loc="left", fontsize=9.5, weight="bold")
        ax.set_xlabel("Surveyed beetle damage, July 2025 (proportion)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(-0.05, 1.0)
        ax.set_ylim(-0.05, 1.0)
    axes[0].set_ylabel("Satellite: share of count area flagged")
    axes[0].legend(frameon=False, fontsize=8, title="GWCT region", title_fontsize=8, loc="upper left")
    fig.suptitle("Satellite vs survey. Grey band = within 0.2 of the survey; hollow = too few clear images. "
                 "Nothing fitted to the survey.", x=0.01, ha="left", fontsize=11, weight="bold")
    plt.tight_layout()
    fig.savefig(out / "fig_beetle_vs_survey.png", dpi=170, facecolor="white")
    log(f"Wrote {out}")


if __name__ == "__main__":
    main()
