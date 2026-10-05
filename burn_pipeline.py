#!/usr/bin/env python3
"""
Heather management detection (burning + cutting) on red grouse count areas.
Sentinel-2 via the Google Earth Engine Python API.

v6, October 2026. Changes from v5 follow Phil's (PW) expert review:

  1. Post-season window now starts 16 April. The burn/cut season closes on
     15 April and early April is often the most intensive period. v5 started
     on 1 April, so the composite straddled the season and early-April scars
     were diluted (they showed up as the "blue" class in v5).
  2. Minimum patch size 0.1 ha (10 pixels). v5 used 0.5 ha, which hid the
     small cuts (sometimes 10 m x 20 m) practised at sites like Crossgill and
     Wemmergill. 0.1 ha is about the floor Sentinel-2 can separate from noise.
  3. Burn/cut timing split removed. Both are practised Oct to 15 Apr, so
     timing cannot separate them. Output is one class: "managed".
  4. Every season is reported, including 2024/25 (heather beetle outbreak)
     and 2025/26. Per-season figures are treated as valid: management is
     genuinely sporadic (intensive one year, nothing the next).
  5. Count area polygons come from a file (definitive boundaries) instead of
     being reconstructed from polylines. Build it with build_count_areas.py.
     Rows with Count_Type "Strip" (ground between Eggleston count areas) are
     measured but kept out of every total; they feed the boundary sensitivity
     summary only.

Usage
  python burn_pipeline.py --polygons data/count_areas_v6.gpkg --diagnose
  python burn_pipeline.py --polygons data/count_areas_v6.gpkg
  python burn_pipeline.py --polygons data/count_areas_v6.gpkg --images --only Edmundbyers
  python burn_pipeline.py --polygons data/count_areas_v6.gpkg --seasons 2024 --exclude-s2c --tag no_s2c
  python burn_pipeline.py --polygons data/count_areas_v6.gpkg --charts-only outputs/results_by_season.csv

Auth (once):  earthengine authenticate
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
from shapely.geometry import box, mapping

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# ---------------------------------------------------------------- settings
DEFAULT_PROJECT = os.environ.get("GEE_PROJECT", "burn-502316")
DEFAULT_SEASONS = list(range(2017, 2026))   # pre-years: 2017 = 2017/18 ... 2025 = 2025/26

BASELINE = ("07-01", "10-01")   # Jul 1 to Sep 30 of the pre-year (end exclusive)
POST     = ("04-16", "07-01")   # Apr 16 to Jun 30 of the post-year (end exclusive)

BANDS = ["B2", "B3", "B4", "B8", "B11", "B12"]
NDVI_MIN = 0.30                 # pixel must have had vegetation to lose
WORLDCOVER_MOORLAND = [20, 30, 90, 100]   # shrubland, grassland, herbaceous wetland, moss/lichen

# Region lookup for charts only. Unknown moors fall into "Other".
REGION = {
    **{m: "North York Moors" for m in ["Danby", "Rosedale", "Westerdale", "Snilesworth"]},
    **{m: "North Pennines" for m in ["Eggleston", "Wemmergill", "Moor House", "Bollihope",
                                     "Crossgill", "Edmundbyers", "Geltsdale"]},
    **{m: "Dales / Nidderdale" for m in ["Grinton", "Swinton", "Ramsgill", "Dallowgill",
                                         "Stags Fell"]},
}
REGION_COLOUR = {"North York Moors": "#c1272d", "North Pennines": "#2b6cb0",
                 "Dales / Nidderdale": "#7a7a7a", "Other": "#b8b8b8"}

# Optional fields carried from the polygon file into the results, with defaults
EXTRA_FIELDS = {"Count_Type": "Core", "GWCT_Region": None}
STRIP = "Strip"                 # Count_Type for ground between count areas: never in totals


def season_name(y):
    return f"{y}/{(y + 1) % 100:02d}"


def log(msg):
    print(msg, flush=True)


# ---------------------------------------------------------------- polygons
def _guess(columns, candidates):
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand in lower:
            return lower[cand]
    return None


def load_polygons(path, moor_field=None, count_field=None):
    """Read count areas, check them, dissolve to one feature per Moor+Count.
    Returns a GeoDataFrame in EPSG:4326 with Moor, Count, area_ha."""
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        sys.exit(f"{path} has no CRS. For a shapefile, the .prj file must sit next to the .shp.")

    gtypes = set(gdf.geom_type.dropna())
    if not gtypes <= {"Polygon", "MultiPolygon"}:
        sys.exit(f"Expected polygons, found {sorted(gtypes)}. "
                 "The 2026 file was polylines and had to be reconstructed; "
                 "the definitive boundaries should be polygons. Check with the supplier.")

    cols = [c for c in gdf.columns if c != "geometry"]
    moor_field = moor_field or _guess(cols, ["moor", "estate", "site", "moor_name"])
    count_field = count_field or _guess(cols, ["count", "count_area", "countarea", "count_name",
                                               "area_name", "name"])
    if not moor_field or not count_field:
        sys.exit(f"Could not identify the moor/count fields. Columns are: {cols}. "
                 "Pass --moor-field and --count-field.")
    log(f"Polygons: {len(gdf)} features | moor field '{moor_field}' | count field '{count_field}' "
        f"| CRS {gdf.crs.to_string()}")

    gdf = gdf.rename(columns={moor_field: "Moor", count_field: "Count"})
    for field, default in EXTRA_FIELDS.items():
        if field not in gdf.columns:
            gdf[field] = default
    gdf = gdf[["Moor", "Count", *EXTRA_FIELDS, "geometry"]]
    gdf["Moor"] = gdf["Moor"].astype(str).str.strip()
    gdf["Count"] = gdf["Count"].astype(str).str.strip()

    gdf = gdf.to_crs(27700)
    gdf["geometry"] = gdf.geometry.buffer(0)                     # repair invalid rings
    n_before = len(gdf)
    gdf = gdf.dissolve(by=["Moor", "Count"], as_index=False)
    if len(gdf) < n_before:
        log(f"  dissolved {n_before} features into {len(gdf)} count areas (multi-part areas merged)")
    gdf["geometry"] = gdf.geometry.simplify(1.0)                  # 1 m, keeps EE payload small
    gdf["area_ha"] = gdf.geometry.area / 1e4

    gdf = gdf.sort_values(["Moor", "Count"]).reset_index(drop=True)
    gdf["uid"] = np.arange(len(gdf))
    core = gdf[gdf.Count_Type != STRIP]
    log(f"  {len(core)} count areas, {core.area_ha.sum():,.0f} ha total, {core.Moor.nunique()} moors"
        + (f" (+ {len(gdf) - len(core)} strip rows, {gdf.area_ha.sum() - core.area_ha.sum():,.0f} ha, "
           "kept out of totals)" if len(core) < len(gdf) else ""))
    return gdf.to_crs(4326)


def select_areas(gdf, only):
    """Count areas whose moor or count name contains any of the comma-separated substrings."""
    if not only:
        return gdf
    keys = [k.strip().lower() for k in only.split(",")]
    mask = gdf.apply(lambda r: any(k in r.Moor.lower() or k in r.Count.lower() for k in keys), axis=1)
    if not mask.any():
        sys.exit(f"--only '{only}' matched no count areas.")
    log(f"  --only kept {mask.sum()} count areas (the seasonal correction still uses all of them)")
    return gdf[mask]


def thumb_region(geom_4326, buffer_m=200):
    """Buffered bounding box, as a coordinate list, for getThumbURL."""
    g = gpd.GeoSeries([geom_4326], crs=4326).to_crs(27700).buffer(buffer_m)
    b = box(*g.total_bounds)
    return json.loads(json.dumps(mapping(gpd.GeoSeries([b], crs=27700).to_crs(4326).iloc[0])))["coordinates"]


# ---------------------------------------------------------------- earth engine
def ee_init(project):
    import ee
    try:
        ee.Initialize(project=project)
    except Exception as e:
        sys.exit(f"Earth Engine init failed for project '{project}': {e}\n"
                 "Run `earthengine authenticate` first, and check the project ID "
                 "(set GEE_PROJECT or pass --project).")
    log(f"Earth Engine initialised (project {project})")
    return ee


def to_fc(ee, gdf):
    feats = []
    for r in gdf.itertuples():
        geo = json.loads(json.dumps(mapping(r.geometry)))
        feats.append(ee.Feature(ee.Geometry(geo, None, False),
                                {"Moor": r.Moor, "Count": r.Count, "uid": int(r.uid)}))
    return ee.FeatureCollection(feats)


class Engine:
    """All Earth Engine logic in one place."""

    def __init__(self, ee, fc, threshold, min_patch_ha, cloud_max, exclude_s2c, aoi_fc=None,
                 post=POST):
        self.ee = ee
        self.fc = fc                      # count areas to report
        # Region for compositing and the seasonal correction: all count areas, so a --only
        # subset gives the same numbers as the full run. The areas themselves, NOT .bounds()
        self.aoi = (fc if aoi_fc is None else aoi_fc).geometry()
        self.threshold = threshold
        self.min_px = max(1, int(round(min_patch_ha * 100)))   # 10 m pixel = 0.01 ha
        self.cloud_max = cloud_max
        self.exclude_s2c = exclude_s2c
        self.post = post                  # after-season window (MM-DD, MM-DD end exclusive)
        wc = ee.Image("ESA/WorldCover/v200/2021").select("Map")
        moor = wc.eq(WORLDCOVER_MOORLAND[0])
        for c in WORLDCOVER_MOORLAND[1:]:
            moor = moor.Or(wc.eq(c))
        self.moorland = moor.rename("moorland")

    def s2(self, start, end, bands=BANDS):
        ee = self.ee
        col = (ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
               .filterBounds(self.aoi).filterDate(start, end))
        if self.exclude_s2c:
            col = col.filter(ee.Filter.neq("SPACECRAFT_NAME", "Sentinel-2C"))
        col = col.select(list(bands) + ["SCL"])  # fixed band order: the S2 archive varies
        cld = (ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
               .filterBounds(self.aoi).filterDate(start, end))
        joined = ee.ImageCollection(ee.Join.saveFirst("cloud_mask").apply(
            primary=col, secondary=cld,
            condition=ee.Filter.equals(leftField="system:index", rightField="system:index")))
        cloud_max = self.cloud_max

        def prep(img):
            prob = ee.Image(img.get("cloud_mask")).select("probability")
            scl = img.select("SCL")
            bad = scl.eq(3).Or(scl.eq(8)).Or(scl.eq(9)).Or(scl.eq(10)).Or(scl.eq(11))
            out = (img.select(list(bands)).divide(10000)
                   .updateMask(prob.lt(cloud_max)).updateMask(bad.Not()))
            out = out.addBands([out.normalizedDifference(["B8", "B12"]).rename("NBR"),
                                out.normalizedDifference(["B8", "B4"]).rename("NDVI")])
            return ee.Image(out.copyProperties(img, ["system:time_start"]))

        return joined.map(prep)

    def windows(self, y):
        return ((f"{y}-{BASELINE[0]}", f"{y}-{BASELINE[1]}"),
                (f"{y + 1}-{self.post[0]}", f"{y + 1}-{self.post[1]}"))

    def deseason(self, dnbr):
        """Subtract the moorland-wide median change for this season (the seasonal
        background: senescence, sun angle, wetness). Only a few % of moorland is
        managed per season, so managed pixels do not move the median."""
        ee = self.ee
        off = dnbr.updateMask(self.moorland).reduceRegion(
            reducer=ee.Reducer.median(), geometry=self.aoi, scale=100,
            maxPixels=1e9, bestEffort=True).get("NBR")
        return dnbr.subtract(ee.Image.constant(ee.Number(off)))

    def season(self, y):
        """Returns dict of images for one season (pre-year y)."""
        ee = self.ee
        (b0, b1), (p0, p1) = self.windows(y)
        base_col = self.s2(b0, b1)
        post_col = self.s2(p0, p1)
        base = base_col.select(["NBR", "NDVI"]).median()
        post_nbr = post_col.select("NBR").median()

        dnbr = self.deseason(base.select("NBR").subtract(post_nbr))
        valid = base.select("NBR").mask().And(post_nbr.mask()).rename("valid")

        raw = (dnbr.gt(self.threshold)
               .And(base.select("NDVI").gt(NDVI_MIN))
               .And(self.moorland)).unmask(0)
        if self.min_px > 1:
            patch = raw.selfMask().toByte().connectedPixelCount(
                maxSize=max(64, self.min_px * 2), eightConnected=True)
            managed = raw.And(patch.gte(self.min_px)).unmask(0)
        else:
            managed = raw
        return {
            "managed": managed.rename("managed"),
            "raw": raw.rename("raw"),
            "valid": valid,
            "dnbr": dnbr.rename("dnbr"),
            "n_clear": post_col.select("NBR").count().rename("n_clear"),
            "base_col": base_col,
            "post_col": post_col,
        }

    def season_stats(self, y):
        ee = self.ee
        s = self.season(y)
        px = ee.Image.pixelArea()
        areas = ee.Image.cat([
            s["managed"].multiply(px).rename("managed_m2"),
            self.moorland.multiply(px).rename("moor_m2"),
            s["valid"].unmask(0).And(self.moorland).multiply(px).rename("visible_m2"),
        ])
        sums = areas.reduceRegions(collection=self.fc, reducer=ee.Reducer.sum(),
                                   scale=10, tileScale=4).getInfo()
        clear = s["n_clear"].unmask(0).reduceRegions(collection=self.fc, reducer=ee.Reducer.mean(),
                                                     scale=10, tileScale=4).getInfo()
        counts = {"baseline_images": s["base_col"].size().getInfo(),
                  "post_images": s["post_col"].size().getInfo()}
        return sums, clear, counts

    def diagnostics(self, y):
        ee = self.ee
        s = self.season(y)
        pct = s["dnbr"].updateMask(self.moorland).reduceRegion(
            reducer=ee.Reducer.percentile([50, 90, 95, 97, 99]),
            geometry=self.aoi, scale=20, maxPixels=1e9, bestEffort=True, tileScale=4).getInfo()
        px = ee.Image.pixelArea()
        steps = ee.Image.cat([
            s["dnbr"].gt(self.threshold).unmask(0).multiply(px).rename("a_above"),
            s["raw"].multiply(px).rename("b_veg_moor"),
            s["managed"].multiply(px).rename("c_patch"),
        ]).reduceRegion(reducer=ee.Reducer.sum(), geometry=self.aoi, scale=10,
                        maxPixels=1e9, tileScale=4).getInfo()
        return {
            "season": season_name(y),
            "baseline_images": s["base_col"].size().getInfo(),
            "post_images": s["post_col"].size().getInfo(),
            **pct,
            "ha_above_threshold": steps["a_above"] / 1e4,
            "ha_after_veg_moorland": steps["b_veg_moor"] / 1e4,
            "ha_after_patch_filter": steps["c_patch"] / 1e4,
        }

    def panel(self, y):
        """False colour post-season composite, managed ground in red, count areas in yellow."""
        ee = self.ee
        s = self.season(y)
        fcol = s["post_col"].median().visualize(bands=["B12", "B8", "B4"], min=0.02, max=0.40)
        red = s["managed"].selfMask().visualize(palette=["ff0000"], opacity=0.6)
        outline = ee.Image().byte().paint(self.fc, 1, 2).visualize(palette=["ffff00"])
        return fcol.blend(red).blend(outline)


# ---------------------------------------------------------------- runs
def run_detection(eng, gdf, seasons, out):
    rows, season_rows = [], []
    meta = gdf.set_index("uid")
    for y in seasons:
        t0 = time.time()
        log(f"Season {season_name(y)} ...")
        for attempt in range(3):
            try:
                sums, clear, counts = eng.season_stats(y)
                break
            except Exception as e:
                if attempt == 2:
                    log(f"  FAILED: {e}")
                    raise
                log(f"  retrying after error: {e}")
                time.sleep(15 * (attempt + 1))
        clear_by_uid = {f["properties"]["uid"]: f["properties"].get("mean") for f in clear["features"]}
        for f in sums["features"]:
            p = f["properties"]
            uid = p["uid"]
            area_ha = meta.loc[uid, "area_ha"]
            managed_ha = (p.get("managed_m2") or 0) / 1e4
            moor_ha = (p.get("moor_m2") or 0) / 1e4
            visible_ha = (p.get("visible_m2") or 0) / 1e4
            rows.append({
                "Moor": p["Moor"], "Count": p["Count"],
                "Count_Type": meta.loc[uid, "Count_Type"], "GWCT_Region": meta.loc[uid, "GWCT_Region"],
                "season": season_name(y), "pre_year": y,
                "area_ha": round(area_ha, 2),
                "moorland_ha": round(moor_ha, 2),
                "visible_ha": round(visible_ha, 2),
                "pct_visible": round(100 * visible_ha / moor_ha, 1) if moor_ha else None,
                "managed_ha": round(managed_ha, 3),
                "pct_of_area": round(100 * managed_ha / area_ha, 3) if area_ha else None,
                "pct_of_moorland": round(100 * managed_ha / moor_ha, 3) if moor_ha else None,
                "post_images_per_pixel": round(clear_by_uid.get(uid) or 0, 1),
            })
        season_rows.append({"season": season_name(y), **counts})
        log(f"  {sum(r['managed_ha'] for r in rows if r['pre_year'] == y and r['Count_Type'] != STRIP):.1f}"
            f" ha managed | "
            f"{counts['post_images']} post-season scenes | {time.time() - t0:.0f}s")
    df = pd.DataFrame(rows)
    return df, pd.DataFrame(season_rows)


def run_diagnostics(eng, seasons, out):
    rows = []
    for y in seasons:
        log(f"Diagnostics {season_name(y)} ...")
        d = eng.diagnostics(y)
        rows.append(d)
        log("  " + " | ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}"
                              for k, v in d.items() if k != "season"))
    df = pd.DataFrame(rows)
    df.to_csv(out / "diagnostics.csv", index=False)
    p95 = [c for c in df.columns if c.endswith("p95")]
    if p95:
        log(f"\nCurrent threshold {eng.threshold}. Median p95 across seasons: "
            f"{df[p95[0]].median():.3f} (range {df[p95[0]].min():.3f} to {df[p95[0]].max():.3f}).")
    log(f"Wrote {out / 'diagnostics.csv'}")
    return df


def fetch_thumb(image, region, dest, label, dims=512):
    """Download one thumbnail, retrying with backoff on HTTP 429 / 5xx. Callers go sequentially."""
    import requests
    for attempt in range(6):
        try:
            url = image.getThumbURL({"region": region, "dimensions": dims, "format": "png"})
            resp = requests.get(url, timeout=300)
            if resp.status_code == 429 or resp.status_code >= 500:
                raise RuntimeError(f"HTTP {resp.status_code}")
            resp.raise_for_status()
            dest.write_bytes(resp.content)
            break
        except Exception as e:
            wait = min(300, 10 * 2 ** attempt)
            log(f"  {label}: {e}; waiting {wait}s")
            time.sleep(wait)
    time.sleep(1)          # sequential and gentle: avoids the 429s the browser hit


def run_images(eng, gdf, seasons, out, results=None, dims=512):
    img_dir = out / "images"
    img_dir.mkdir(parents=True, exist_ok=True)
    panels = {y: eng.panel(y) for y in seasons}
    gdf = gdf[gdf.Count_Type != STRIP]
    total = len(gdf) * len(seasons)
    n = 0
    for r in gdf.itertuples():
        tag = f"{r.Moor}_{r.Count}".replace(" ", "")
        region = thumb_region(r.geometry)
        (img_dir / tag).mkdir(exist_ok=True)
        for y in seasons:
            n += 1
            dest = img_dir / tag / f"{y}_{y + 1}.png"
            if dest.exists():
                continue
            fetch_thumb(panels[y], region, dest, f"{tag} {season_name(y)}", dims)
            log(f"  [{n}/{total}] {tag} {season_name(y)}")
        if results is not None:
            contact_sheet(img_dir / tag, r.Moor, r.Count, r.area_ha, seasons, results, out)
    log(f"Images in {img_dir}")


# ---------------------------------------------------------------- charts
def _style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlelocation": "left",
                         "axes.titleweight": "bold"})


def contact_sheet(folder, moor, count, area_ha, seasons, results, out):
    import matplotlib.image as mpimg
    files = [(y, folder / f"{y}_{y + 1}.png") for y in seasons]
    files = [(y, f) for y, f in files if f.exists()]
    if not files:
        return
    _style()
    ncol = 5
    nrow = int(np.ceil(len(files) / ncol))
    first = mpimg.imread(files[0][1])
    aspect = first.shape[0] / first.shape[1]                     # rows follow the image shape
    head = 0.8                                                    # inches for title and key
    fig_h = nrow * (3.2 * aspect + 0.75) + head
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, fig_h))
    axes = np.atleast_1d(axes).ravel()
    sub = results[(results.Moor == moor) & (results.Count == count)].set_index("pre_year")
    for ax in axes:
        ax.axis("off")
    for ax, (y, f) in zip(axes, files):
        ax.imshow(mpimg.imread(f))
        ax.set_title(season_name(y), fontsize=11)
        if y in sub.index:
            ha, pct = sub.loc[y, "managed_ha"], sub.loc[y, "pct_of_area"]
            ax.text(0.5, -0.04, f"{ha:.1f} ha ({pct:.1f}%)", transform=ax.transAxes, va="top",
                    ha="center", fontsize=9.5, color="#c1272d" if pct >= 10 else "#333")
    fig.suptitle(f"{moor}, {count} ({area_ha:.0f} ha): managed ground by season",
                 x=0.01, y=1 - 0.12 / fig_h, ha="left", va="top", fontsize=13, weight="bold")
    fig.text(0.01, 1 - 0.48 / fig_h, "Sentinel-2 false colour (SWIR, NIR, red). "
             "Red = detected burning or cutting. Yellow = count area.", fontsize=9, color="#555", va="top")
    plt.tight_layout(rect=[0, 0, 1, 1 - head / fig_h])
    dest = out / "figures" / "contact_sheets"
    dest.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest / f"{folder.name}.png", dpi=160, facecolor="white")
    plt.close(fig)


def boundary_sensitivity(df, out):
    """Where a supplied polygon was split into count areas plus a strip (Eggleston), compare
    % managed for the count areas, the strip, and the whole supplied polygon, by moor and region."""
    moors = df.loc[df.Count_Type == STRIP, "Moor"].unique()
    if not len(moors):
        return
    rows = []
    for scope, d in [*((m, df[df.Moor == m]) for m in moors),
                     *((r, df[df.Region == r]) for r in df.loc[df.Moor.isin(moors), "Region"].unique())]:
        for label, sub in [("count areas", d[d.Count_Type != STRIP]),
                           ("strip between counts", d[d.Count_Type == STRIP]),
                           ("count areas + strip", d)]:
            s = sub.groupby("season").agg(m=("managed_ha", "sum"), a=("area_ha", "sum"))
            rows.append({"scope": scope, "boundary": label, "area_ha": round(s.a.iloc[0], 1),
                         "mean_annual_pct": round(100 * s.m.sum() / s.a.sum(), 2),
                         **(100 * s.m / s.a).round(2).to_dict()})
    pd.DataFrame(rows).to_csv(out / "summary_boundary_sensitivity.csv", index=False)
    log(f"Boundary sensitivity in {out / 'summary_boundary_sensitivity.csv'}")


def make_charts(df, gdf, out, figs=None, titles=True):
    """titles=False and figs=<folder> give caption-ready versions for the report."""
    _style()
    figs = Path(figs) if figs else out / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    df = df.copy()
    df["Region"] = df.Moor.map(REGION).fillna("Other")
    if "Count_Type" in df.columns:
        boundary_sensitivity(df, out)
        df = df[df.Count_Type != STRIP]
    gdf = gdf[gdf.Count_Type != STRIP]
    seasons = sorted(df.season.unique(), key=lambda s: int(s[:4]))

    # 1. Total managed area per season
    tot = df.groupby("season").agg(managed_ha=("managed_ha", "sum"), area_ha=("area_ha", "sum"))
    tot = tot.loc[seasons]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    x = np.arange(len(seasons))
    ax.bar(x, tot.managed_ha, color="#c1272d", zorder=3)
    for i, (ha, a) in enumerate(zip(tot.managed_ha, tot.area_ha)):
        ax.text(i, ha + tot.managed_ha.max() * 0.015, f"{ha:.0f} ha\n{100 * ha / a:.1f}%",
                ha="center", fontsize=8)
    if "2021/22" in seasons:
        i = seasons.index("2021/22") - 0.5
        ax.axvline(i, color="#2b6cb0", ls="--", lw=1.3)
        ax.text(i + 0.05, tot.managed_ha.max() * 1.12, "Deep-peat burning\nregulation, May 2021",
                color="#2b6cb0", fontsize=8.5, va="top")
    if "2024/25" in seasons:
        j = seasons.index("2024/25")
        ax.annotate("heather beetle\noutbreak", xy=(j, tot.managed_ha.iloc[j]),
                    xytext=(j - 1.7, tot.managed_ha.max() * 0.85), fontsize=8.5, color="#555",
                    arrowprops=dict(arrowstyle="->", color="#555"))
    ax.set_xticks(x)
    ax.set_xticklabels(seasons)
    ax.set_ylim(0, tot.managed_ha.max() * 1.22)
    ax.set_ylabel("Managed area, burnt or cut (ha)")
    if titles:
        ax.set_title("Burnt or cut ground within count areas, by season")
    ax.grid(axis="y", alpha=0.25, zorder=0)
    plt.tight_layout()
    fig.savefig(figs / "fig_season_totals.png", dpi=200)
    plt.close(fig)

    # 2. Count area x season heatmap (shows the sporadic, intensive pattern)
    piv = df.pivot_table(index=["Moor", "Count"], columns="season", values="pct_of_area").reindex(columns=seasons)
    piv = piv.loc[piv.mean(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(0.62 * len(seasons) + 3.0, 0.24 * len(piv) + 1.3))   # prints legibly at A4 width
    vmax = max(10, np.nanpercentile(piv.values, 98))
    im = ax.imshow(piv.values, cmap="YlOrRd", vmin=0, vmax=vmax, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v) and v >= 0.5:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                        color="white" if v > vmax * 0.6 else "#222")
    ax.set_xticks(range(len(seasons)))
    ax.set_xticklabels(seasons, fontsize=8.5, rotation=45, ha="right")
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([f"{m}, {c}" for m, c in piv.index], fontsize=8.5)
    ax.spines[:].set_visible(False)
    plt.colorbar(im, ax=ax, shrink=0.5, pad=0.02, label="% of count area burnt or cut")
    if titles:
        ax.set_title("Managed ground by count area and season (%)")
    plt.tight_layout()
    fig.savefig(figs / "fig_heatmap.png", dpi=200)
    plt.close(fig)

    # 3. Mean annual % by moor
    bym = df.groupby(["Moor", "Region"]).agg(m=("managed_ha", "sum"), a=("area_ha", "sum")).reset_index()
    bym["pct"] = 100 * bym.m / bym.a
    bym = bym.sort_values("pct")
    fig, ax = plt.subplots(figsize=(8, 0.32 * len(bym) + 1.4))
    ax.barh(bym.Moor, bym.pct, color=bym.Region.map(REGION_COLOUR), zorder=3)
    for i, v in enumerate(bym.pct):
        ax.text(v + bym.pct.max() * 0.01, i, f"{v:.1f}%", va="center", fontsize=8)
    ax.set_xlabel("Mean annual % of count area burnt or cut")
    if titles:
        ax.set_title("Management intensity by moor")
    from matplotlib.patches import Patch
    present = [r for r in REGION_COLOUR if r in set(bym.Region)]
    ax.legend(handles=[Patch(color=REGION_COLOUR[r], label=r) for r in present],
              frameon=False, fontsize=8.5, loc="lower right")
    ax.grid(axis="x", alpha=0.25, zorder=0)
    plt.tight_layout()
    fig.savefig(figs / "fig_by_moor.png", dpi=200)
    plt.close(fig)

    # 4. Map: one labelled, proportional circle per moor (answers PW: "can't see the proportion")
    g = gdf.to_crs(27700).copy()
    cent = g.dissolve(by="Moor").centroid
    pct = bym.set_index("Moor").pct
    region_of = bym.set_index("Moor").Region
    fig, ax = plt.subplots(figsize=(8.5, 8))
    g.plot(ax=ax, facecolor="#eeeeee", edgecolor="#999", lw=0.5)
    pad = 8000
    xmin, ymin, xmax, ymax = g.total_bounds
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad * 2)                       # headroom for region names
    ax.set_aspect("equal")
    ax.axis("off")
    vmax = max(pct.max(), 1)
    circles = []
    for moor, pt in cent.items():
        v = pct.get(moor, 0)
        size = 80 + 900 * v / vmax
        ax.scatter(pt.x, pt.y, s=size, color=plt.cm.YlOrRd(0.15 + 0.85 * v / vmax),
                   edgecolor="#333", lw=0.6, zorder=3, alpha=0.9)
        circles.append((moor, pt, v, np.sqrt(size / np.pi)))          # radius in points
    # 10 km scale bar, bottom left
    sx, sy = xmin - pad * 0.6, ymin - pad * 0.6
    ax.plot([sx, sx + 10000], [sy, sy], color="#333", lw=2)
    ax.text(sx + 5000, sy + 1200, "10 km", ha="center", fontsize=8, color="#333")

    # Labels: try positions round each circle, keep the first that hits no other label or circle
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    px = fig.dpi / 72
    centres = {m: (ax.transData.transform((p.x, p.y)), r * px) for m, p, _, r in circles}

    def clashes(bb, own, placed):
        if any(bb.overlaps(o) for o in placed):
            return True
        for m, (c, rad) in centres.items():
            if m == own:
                continue
            nx, ny = min(max(c[0], bb.x0), bb.x1), min(max(c[1], bb.y0), bb.y1)
            if np.hypot(c[0] - nx, c[1] - ny) < rad:
                return True
        return False

    spots = [(1, 0, "left", "center"), (-1, 0, "right", "center"), (0, 1, "center", "bottom"),
             (0, -1, "center", "top"), (0.7, 0.7, "left", "bottom"), (-0.7, 0.7, "right", "bottom"),
             (0.7, -0.7, "left", "top"), (-0.7, -0.7, "right", "top")]
    placed = []
    # Region names first, above each cluster, so moor labels avoid them
    for region, moors in region_of.groupby(region_of):
        pts = [cent[m] for m in moors.index]
        x, y = np.mean([p.x for p in pts]), max(p.y for p in pts)
        t = ax.annotate(region, (x, y), xytext=(0, 34), textcoords="offset points", ha="center",
                        fontsize=10, color="#888", style="italic", zorder=2)
        placed.append(t.get_window_extent(rend).expanded(1.05, 1.2))
    for moor, pt, v, r in sorted(circles, key=lambda c: -c[3]):
        for i, (dx, dy, ha, va) in enumerate(spots):
            t = ax.annotate(f"{moor}\n{v:.1f}%/yr", (pt.x, pt.y), xytext=((r + 3) * dx, (r + 3) * dy),
                            textcoords="offset points", ha=ha, va=va, fontsize=8, zorder=4,
                            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1))
            bb = t.get_window_extent(rend).expanded(1.04, 1.08)
            if not clashes(bb, moor, placed) or i == len(spots) - 1:
                break
            t.remove()
        placed.append(bb)
    if titles:
        ax.set_title("Mean annual % of count area burnt or cut, by moor\n(circle size and colour = %)")
    plt.tight_layout()
    fig.savefig(figs / "fig_map.png", dpi=200, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)

    # 5. Summary tables
    by_region = df.groupby("Region").agg(m=("managed_ha", "sum"), a=("area_ha", "sum"))
    by_region["mean_annual_pct"] = (100 * by_region.m / by_region.a).round(2)
    by_region[["mean_annual_pct"]].to_csv(out / "summary_by_region.csv")
    bym[["Moor", "Region", "pct"]].rename(columns={"pct": "mean_annual_pct"}).round(2) \
        .sort_values("mean_annual_pct", ascending=False).to_csv(out / "summary_by_moor.csv", index=False)
    tot.assign(pct=(100 * tot.managed_ha / tot.area_ha).round(2)).to_csv(out / "summary_by_season.csv")
    log(f"Charts in {figs}; summaries in {out}")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--polygons", required=True, help="count area polygons (shp, gpkg, geojson, kml)")
    ap.add_argument("--moor-field")
    ap.add_argument("--count-field")
    ap.add_argument("--only", help="comma-separated substrings to restrict count areas, e.g. 'Edmundbyers,Moor House'")
    ap.add_argument("--project", default=DEFAULT_PROJECT)
    ap.add_argument("--seasons", default=",".join(map(str, DEFAULT_SEASONS)),
                    help="pre-years, e.g. 2018,2024 or 2017-2025")
    ap.add_argument("--threshold", type=float, default=0.17, help="deseasoned dNBR threshold (v5 p95)")
    ap.add_argument("--min-patch-ha", type=float, default=0.1)
    ap.add_argument("--cloud-max", type=int, default=40)
    ap.add_argument("--post-window", default=",".join(POST),
                    help="after-season window MM-DD,MM-DD (end exclusive); v5 used 04-01,06-01")
    ap.add_argument("--exclude-s2c", action="store_true", help="drop Sentinel-2C scenes (2024/25 sensitivity test)")
    ap.add_argument("--diagnose", action="store_true", help="threshold percentiles and pixel-loss per step only")
    ap.add_argument("--images", action="store_true", help="also download a panel per count area per season")
    ap.add_argument("--charts-only", help="skip Earth Engine; rebuild charts from this results CSV")
    ap.add_argument("--out", default="outputs")
    ap.add_argument("--tag", default="", help="suffix for the output folder, e.g. no_s2c")
    a = ap.parse_args()

    if "-" in a.seasons:
        lo, hi = map(int, a.seasons.split("-"))
        seasons = list(range(lo, hi + 1))
    else:
        seasons = [int(s) for s in a.seasons.split(",")]

    out = Path(a.out + (f"_{a.tag}" if a.tag else ""))
    out.mkdir(parents=True, exist_ok=True)
    gdf_all = load_polygons(a.polygons, a.moor_field, a.count_field)
    gdf = select_areas(gdf_all, a.only)
    gdf.to_file(out / "count_areas_used.geojson", driver="GeoJSON")

    if a.charts_only:
        make_charts(pd.read_csv(a.charts_only), gdf, out)
        return

    ee = ee_init(a.project)
    eng = Engine(ee, to_fc(ee, gdf), a.threshold, a.min_patch_ha, a.cloud_max, a.exclude_s2c,
                 aoi_fc=None if gdf is gdf_all else to_fc(ee, gdf_all),
                 post=tuple(a.post_window.split(",")))
    config = {k: v for k, v in vars(a).items()}
    config.update(seasons=seasons, baseline_window=BASELINE, post_window=eng.post,
                  ndvi_min=NDVI_MIN, worldcover_classes=WORLDCOVER_MOORLAND)
    (out / "run_config.json").write_text(json.dumps(config, indent=2))

    if a.diagnose:
        run_diagnostics(eng, seasons, out)
        return

    df, sdf = run_detection(eng, gdf, seasons, out)
    df.to_csv(out / "results_by_season.csv", index=False)
    sdf.to_csv(out / "scene_counts.csv", index=False)
    log(f"Wrote {out / 'results_by_season.csv'} ({len(df)} rows)")
    make_charts(df, gdf, out)

    if a.images:
        run_images(eng, gdf, seasons, out, results=df)


if __name__ == "__main__":
    main()
