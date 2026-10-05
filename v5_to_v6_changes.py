"""How much did each v5 -> v6 change move the results?

v6 changed three things at once after PW's review: the count area boundaries
(Eleanor's polygons replace v5's reconstruction), the after-season window
(16 Apr to 30 Jun instead of 1 Apr to 31 May) and the minimum patch (0.1 ha
instead of 0.5 ha). This script separates them.

  python v5_to_v6_changes.py --build      # rebuild v5's polygons from the transect lines
  # then run the pipeline once per entry in RUNS (printed by --build; ~15 min in parallel)
  python v5_to_v6_changes.py --summarise  # outputs/v5_to_v6/*.csv and figure

v5's polygons are rebuilt exactly as v5 did it: the closed ring where a line was
drawn as one, otherwise the convex hull of the count's lines. Areas match v5's
to within 0.5%.
"""
import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
from shapely.ops import linemerge

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LINES = "data/source/PW_RGC_2026_lines/PW RGC 2026.shp"
V5_AREAS = "data/count_areas_v5_reconstructed.gpkg"
OUT = Path("outputs/v5_to_v6")

V5_WINDOW, V6_WINDOW = "04-01,06-01", "04-16,07-01"
# tag: (label, boundaries, after window, min patch ha)
RUNS = {
    "abl_v5b_v5s": ("v5 boundaries, v5 settings", "v5", V5_WINDOW, 0.5),
    "abl_v5b_v6s": ("v5 boundaries, v6 settings", "v5", V6_WINDOW, 0.1),
    "abl_v6b_v5s": ("v6 boundaries, v5 settings", "v6", V5_WINDOW, 0.5),
    "abl_v6b_date": ("v6 boundaries, 16 April window only", "v6", V6_WINDOW, 0.5),
    "abl_v6b_patch": ("v6 boundaries, 0.1 ha patch only", "v6", V5_WINDOW, 0.1),
    "abl_v6b_start": ("v6 boundaries, 16 April start only", "v6", "04-16,06-01", 0.5),
    "abl_v6b_start01": ("v6 boundaries, 16 April start and 0.1 ha", "v6", "04-16,06-01", 0.1),
    "": ("v6 boundaries, v6 settings (main run)", "v6", V6_WINDOW, 0.1),
}
UNCHANGED = ["Dallowgill|Bishops", "Ramsgill|Raygill", "Swinton|Potts", "Grinton|Long Gill",
             "Snilesworth|Lodge", "Danby|Low Moor", "Moor House|Behind House"]
V5_SEASONS = [f"{y}/{str(y + 1)[2:]}" for y in range(2017, 2025)]


def build():
    g = gpd.read_file(LINES).dropna(subset=["Moor"])
    rows = []
    for (moor, count), grp in g.groupby(["Moor", "Count"]):
        geom = grp.union_all().convex_hull
        if len(grp) == 1:   # one line drawn as a (near-)closed ring: keep its shape
            ln = grp.geometry.iloc[0]
            ln = linemerge(ln) if ln.geom_type == "MultiLineString" else ln
            if ln.geom_type == "LineString" and len(ln.coords) > 3:
                ring = Polygon(ln.coords).buffer(0)
                if ring.area > 0.9 * geom.area:
                    geom = ring
        rows.append({"Moor": moor, "Count": count, "Count_Type": "Core", "geometry": geom})
    out = gpd.GeoDataFrame(rows, crs=g.crs)
    out.to_file(V5_AREAS, driver="GPKG")
    v5 = pd.read_csv("reference/results_v5.csv").groupby(["Moor", "Count"]).area_ha.first()
    chk = out.set_index(["Moor", "Count"]).area / 1e4
    print(f"Wrote {V5_AREAS}: {len(out)} count areas, {chk.sum():,.0f} ha (v5: {v5.sum():,.0f} ha); "
          f"largest area difference {(chk / v5 - 1).abs().max():.1%}")
    for tag, (label, b, win, patch) in RUNS.items():
        if tag:
            poly = V5_AREAS if b == "v5" else "data/count_areas_v6.gpkg"
            print(f"python burn_pipeline.py --polygons {poly} --post-window {win} "
                  f"--min-patch-ha {patch} --tag {tag}")


def load(tag):
    df = pd.read_csv(Path("outputs" + (f"_{tag}" if tag else "")) / "results_by_season.csv")
    if "Count_Type" in df:
        df = df[df.Count_Type != "Strip"]
    df["key"] = df.Moor + "|" + df.Count
    return df


def summarise():
    OUT.mkdir(parents=True, exist_ok=True)
    runs = {tag: load(tag) for tag in RUNS}
    common = set.intersection(*[set(d.key) for d in runs.values()])   # 27: Glaisdale is v6 only
    v5 = pd.read_csv("reference/results_v5.csv")
    v5["key"] = v5.Moor + "|" + v5.Count

    rows = []
    for tag, (label, b, win, patch) in RUNS.items():
        d = runs[tag]
        d8 = d[d.key.isin(common) & d.season.isin(V5_SEASONS)]
        un = d8[d8.key.isin(UNCHANGED)]
        seas = d8.groupby("season").managed_ha.sum()
        rows.append({
            "run": label, "boundaries": b, "after_window": win, "min_patch_ha": patch,
            "area_ha": d8.groupby("key").area_ha.first().sum(),
            "mean_ha_per_season": seas.mean(),
            "pct_per_year": 100 * d8.managed_ha.sum() / d8.area_ha.sum(),
            "pct_seasons_zero": 100 * (d8.managed_ha == 0).mean(),
            "unchanged7_total_ha": un.managed_ha.sum(),
            "pre2021_ha": seas[V5_SEASONS[:4]].mean(), "post2021_ha": seas[V5_SEASONS[4:7]].mean(),
            "s2024_ha": seas[V5_SEASONS[7]],
        })
    v58 = v5[v5.key.isin(common)]
    vs = v58.groupby("season")
    rows.insert(0, {"run": "v5 as reported (burn class only)", "boundaries": "v5",
                    "after_window": V5_WINDOW, "min_patch_ha": 0.5,
                    "area_ha": v58.groupby("key").area_ha.first().sum(),
                    "mean_ha_per_season": vs.burn_ha.sum().mean(),
                    "pct_per_year": 100 * v58.burn_ha.sum() / v58.area_ha.sum(),
                    "pct_seasons_zero": 100 * (v58.burn_ha == 0).mean(),
                    "unchanged7_total_ha": v58[v58.key.isin(UNCHANGED)].burn_ha.sum(),
                    "pre2021_ha": vs.burn_ha.sum()[V5_SEASONS[:4]].mean(),
                    "post2021_ha": vs.burn_ha.sum()[V5_SEASONS[4:7]].mean(),
                    "s2024_ha": vs.burn_ha.sum()[V5_SEASONS[7]]})
    tab = pd.DataFrame(rows)
    tab.round(2).to_csv(OUT / "v5_to_v6_summary.csv", index=False)
    print(tab.drop(columns=["after_window"]).round(1).to_string(index=False))

    # Per count area: how much does each change move the mean annual %?
    def pct(tag, keys=common):
        d = runs[tag]
        d = d[d.key.isin(keys) & d.season.isin(V5_SEASONS)]
        return d.groupby("key").apply(lambda x: 100 * x.managed_ha.sum() / x.area_ha.sum())

    per = pd.DataFrame({
        "v5_settings_v5_boundaries": pct("abl_v5b_v5s"),
        "v5_settings_v6_boundaries": pct("abl_v6b_v5s"),
        "date_only": pct("abl_v6b_date"),
        "patch_only": pct("abl_v6b_patch"),
        "v6": pct(""),
    }).round(2)
    per.to_csv(OUT / "v5_to_v6_by_count_area.csv")

    # Seasons: the date change by season (v6 boundaries, 0.5 ha patch)
    def by_season(tag):
        d = runs[tag]
        return d[d.key.isin(common)].groupby("season").managed_ha.sum()

    seas = pd.DataFrame({"1 Apr window": by_season("abl_v6b_v5s"),
                         "16 Apr window": by_season("abl_v6b_date"),
                         "16 Apr + 0.1 ha (v6)": by_season("")}).round(1)
    seas.to_csv(OUT / "v5_to_v6_by_season.csv")
    print(seas.to_string())

    # Figure: step chart from v5 settings to v6, on v6 boundaries (27 common count areas)
    t = tab.set_index("run")
    steps = [("v5 settings\nv5 boundaries", t.loc["v5 boundaries, v5 settings", "mean_ha_per_season"] /
              t.loc["v5 boundaries, v5 settings", "area_ha"] * 100),
             ("v5 settings\nv6 boundaries", t.loc["v6 boundaries, v5 settings", "pct_per_year"]),
             ("+ 16 April\nwindow", t.loc["v6 boundaries, 16 April window only", "pct_per_year"]),
             ("+ 0.1 ha\npatch (v6)", t.loc["v6 boundaries, v6 settings (main run)", "pct_per_year"])]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    xs = np.arange(len(steps))
    vals = [v for _, v in steps]
    ax.bar(xs, vals, width=0.55, color=["#9aa5b1", "#9aa5b1", "#2b6cb0", "#c1272d"])
    for x, v in zip(xs, vals):
        ax.text(x, v + 0.08, f"{v:.1f}%", ha="center", va="bottom", fontsize=10)
    ax.set_xticks(xs, [s for s, _ in steps], fontsize=9.5)
    ax.set_ylabel("Mean % of count area\nburnt or cut per season", fontsize=9.5)
    ax.set_ylim(0, max(vals) * 1.2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_title("27 count areas common to v5 and v6, seasons 2017/18 to 2024/25", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(OUT / "fig_v5_to_v6.png", dpi=200, bbox_inches="tight")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--summarise", action="store_true")
    a = ap.parse_args()
    if a.build:
        build()
    if a.summarise:
        summarise()
