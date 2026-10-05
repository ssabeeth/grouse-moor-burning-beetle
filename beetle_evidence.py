#!/usr/bin/env python3
"""
Photo evidence for heather beetle damage: show the imagery behind the satellite estimate, at
the count areas where it is furthest from the July 2025 survey and, for contrast, closest.

Sites are chosen by rule from outputs/beetle_check/beetle_vs_survey.csv (run beetle_check.py
first): among count areas with enough clear images, the 5 largest and 5 smallest gaps between
the satellite estimate (share flagged) and the survey.

One row per count area, six Sentinel-2 panels, all the same shape:

  True colour, late summer (Jul to Sep) 2024, 2025 and 2026: what the eye sees
  Browning   change in NBR, summer 2025 vs 2024, beyond the moor-wide change (brown = browner)
  Greying    change in colour saturation, late vs early summer 2025, beyond the moor-wide
             change (grey = greyer)
  Flagged    pixels the estimate counts: browned by more than 0.10 or greyed by more than 0.05
             beyond the moor-wide change (red)

Ground burnt or cut in 2024/25 or 2025/26, plus a 20 m margin, is black on the maps and left
out of the estimate. Count area outlined (yellow on photos, black on maps).

Usage
  python beetle_evidence.py --polygons data/count_areas_v6.gpkg
"""

import argparse
import json
import shutil
from pathlib import Path

import geopandas as gpd
import matplotlib
import pandas as pd
from shapely.geometry import box, mapping

from beetle_check import FLAG, composite
from burn_pipeline import DEFAULT_PROJECT, Engine, ee_init, fetch_thumb, load_polygons, log, to_fc

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

N_SITES = 5
BUFFER_M = 20                                       # margin around burns and cuts, as beetle_check v1
ASPECT, PAD_M = 1.5, 250                            # every panel: width / height, padding round the area
PANELS = ["tc_2024", "tc_2025", "tc_2026", "browning", "greying", "flagged"]
TITLES = {"tc_2024": "True colour, Jul-Sep 2024", "tc_2025": "True colour, Jul-Sep 2025",
          "tc_2026": "True colour, Jul-Sep 2026", "browning": "Browning, summer 2025 vs 2024",
          "greying": "Greying, late vs early summer 2025", "flagged": "Counted as beetle damage"}
TRUE_COLOUR = dict(bands=["B4", "B3", "B2"], min=0.01, max=0.13, gamma=1.3)
BROWN_PALETTE = ["8c510a", "d8b365", "f6e8c3", "f5f5f5", "c7eae5", "5ab4ac", "01665e"]   # brown ... teal
GREY_PALETTE = ["3a3a3a", "8a8a8a", "d4d4d4", "f7f7f7", "d8daeb", "998ec3", "542788"]    # grey ... purple
NOTE = ("Maps: brown = browned beyond the moor-wide change; grey = greyed beyond it; teal / purple = the opposite; "
        "red = counted as beetle damage; light grey = moorland not counted;\n"
        "black = burnt or cut in 2024/25 or 2025/26 plus a 20 m margin (left out); white = not moorland. "
        "Survey: July 2025 proportion damaged. Satellite: share of usable count area flagged red. "
        "Nothing is fitted to the survey.")


def fixed_region(geom_4326):
    """Box round the count area, padded and widened to a fixed aspect, as a coordinate list."""
    x0, y0, x1, y1 = gpd.GeoSeries([geom_4326], crs=4326).to_crs(27700).total_bounds
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    w, h = x1 - x0 + 2 * PAD_M, y1 - y0 + 2 * PAD_M
    w, h = max(w, h * ASPECT), max(h, w / ASPECT)
    b = gpd.GeoSeries([box(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)], crs=27700).to_crs(4326).iloc[0]
    return json.loads(json.dumps(mapping(b)))["coordinates"]


def build_images(ee, eng):
    tc = {y: eng.s2(f"{y}-07-01", f"{y}-10-01").select(["B4", "B3", "B2"]).median() for y in (2024, 2025, 2026)}
    late24, early25, late25 = (composite(ee, eng, *w) for w in [("2024-07-01", "2024-10-01"),
                                                                 ("2025-04-16", "2025-07-01"),
                                                                 ("2025-07-01", "2025-10-01")])
    diff = ee.Image.cat([late25.select("NBR").subtract(late24.select("NBR")),
                         late25.select("Sat").subtract(early25.select("Sat"))]).rename(list(FLAG))
    managed = eng.season(2024)["managed"].Or(eng.season(2025)["managed"])
    drop = managed.focalMax(radius=BUFFER_M, units="meters")
    keep = eng.moorland.And(drop.Not())
    off = diff.updateMask(keep).reduceRegion(reducer=ee.Reducer.median(), geometry=eng.aoi, scale=20,
                                             maxPixels=1e9, bestEffort=True, tileScale=4).getInfo()
    log("Moor-wide median change removed: " + ", ".join(f"{k} {v:+.3f}" for k, v in off.items()))
    dese = diff.subtract(ee.Image.constant([off[k] for k in FLAG]).rename(list(FLAG)))
    nbr, sat = (dese.select(k) for k in FLAG)
    flagged = None
    for k, thr in FLAG.items():
        f = dese.select(k).lt(thr)
        flagged = f if flagged is None else flagged.Or(f)

    white = ee.Image.constant(1).visualize(palette=["ffffff"])
    black = drop.selfMask().visualize(palette=["000000"])
    edge_y = ee.Image().byte().paint(eng.fc, 1, 2).visualize(palette=["ffff00"])
    edge_k = ee.Image().byte().paint(eng.fc, 1, 1).visualize(palette=["000000"])

    def on_map(vis):
        return white.blend(vis).blend(black).blend(edge_k)

    imgs = {f"tc_{y}": im.visualize(**TRUE_COLOUR).blend(edge_y) for y, im in tc.items()}
    imgs["browning"] = on_map(nbr.updateMask(eng.moorland).visualize(min=-0.25, max=0.25, palette=BROWN_PALETTE))
    imgs["greying"] = on_map(sat.updateMask(eng.moorland).visualize(min=-0.15, max=0.15, palette=GREY_PALETTE))
    imgs["flagged"] = on_map(keep.selfMask().visualize(palette=["e3e3e3"])
                             .blend(flagged.updateMask(keep).selfMask().visualize(palette=["d7301f"])))
    return imgs


def pick_sites(comp):
    ok = comp[comp.enough_images].assign(gap=lambda d: (d.satellite_estimate - d.damage).abs())
    return {"disagree": ok.nlargest(N_SITES, "gap"), "agree": ok.nsmallest(N_SITES, "gap")}


def sheet(name, sites, img_dir, out):
    fig, axes = plt.subplots(len(sites), len(PANELS), figsize=(2.9 * len(PANELS), 2.15 * len(sites) + 1.3))
    for i, r in enumerate(sites.itertuples()):
        tag = f"{r.Moor}_{r.Count}".replace(" ", "")
        for j, p in enumerate(PANELS):
            ax = axes[i, j]
            ax.imshow(mpimg.imread(img_dir / tag / f"{p}.png"))
            ax.set_xticks([])
            ax.set_yticks([])
            if i == 0:
                ax.set_title(TITLES[p], fontsize=9, loc="left")
        axes[i, 0].set_ylabel(f"{r.Moor}, {r.Count}\nsurvey {r.damage:.1f} | satellite {r.satellite_estimate:.2f}",
                              fontsize=9, weight="bold")
    title = ("Heather beetle: the 5 count areas where satellite and survey differ most" if name == "disagree"
             else "Heather beetle: the 5 count areas where satellite and survey agree best, for contrast")
    fig.suptitle(title, x=0.01, ha="left", fontsize=12.5, weight="bold")
    fig.text(0.01, 0.008, NOTE, fontsize=7.5, color="#444", va="bottom")
    plt.tight_layout(rect=[0, 0.045, 1, 0.97])
    fig.savefig(out / f"fig_evidence_{name}.png", dpi=140, facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--polygons", default="data/count_areas_v6.gpkg")
    ap.add_argument("--project", default=DEFAULT_PROJECT)
    ap.add_argument("--threshold", type=float, default=0.17, help="for the managed mask")
    ap.add_argument("--out", default="outputs/beetle_check")
    a = ap.parse_args()
    out = Path(a.out)
    img_dir = out / "evidence"
    shutil.rmtree(img_dir, ignore_errors=True)       # panel regions depend on the settings above
    sheets = pick_sites(pd.read_csv(out / "beetle_vs_survey.csv"))
    for name, s in sheets.items():
        log(f"{name}: " + ", ".join(f"{r.Count} ({r.damage:.1f} vs {r.satellite_estimate:.2f})" for r in s.itertuples()))

    gdf = load_polygons(a.polygons)
    ee = ee_init(a.project)
    eng = Engine(ee, to_fc(ee, gdf), a.threshold, 0.1, 40, False)
    imgs = build_images(ee, eng)

    sites = pd.concat(sheets.values())
    for k, r in enumerate(sites.itertuples(), 1):
        geom = gdf[(gdf.Moor == r.Moor) & (gdf.Count == r.Count)].geometry.iloc[0]
        tag = f"{r.Moor}_{r.Count}".replace(" ", "")
        (img_dir / tag).mkdir(parents=True, exist_ok=True)
        region = fixed_region(geom)
        for p in PANELS:
            fetch_thumb(imgs[p], region, img_dir / tag / f"{p}.png", f"{tag} {p}", dims=480)
        log(f"  [{k}/{len(sites)}] {r.Moor}, {r.Count}")
    for name, s in sheets.items():
        sheet(name, s, img_dir, out)
    log(f"Evidence sheets in {out}")


if __name__ == "__main__":
    main()
