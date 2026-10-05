#!/usr/bin/env python3
"""
Build the named count area file for v6 from the two supplied layers.

  Eleanor's definitive polygons (Oct 2026): 22 polygons, correct boundaries, but no names.
  The PW_RGC_2026 lines (Jun 2026): the count transects, with Moor and Count names.

Each polygon takes the name of the transects it contains. Two exceptions:

  Danby, Glaisdale  Its two transects are unnamed in PW_RGC_2026. Confirmed as
                    Glaisdale by location (NZ 73 03) and the 2025 beetle sheet.
  Eggleston         Eleanor drew all 7 count areas as one 1,928 ha block. Each count
                    area is rebuilt with her own rule (the two transects plus 150 m),
                    clipped to her block. What is left (the strips between count areas)
                    is kept as one extra row, Count_Type "Strip". The pipeline leaves it
                    out of all totals; counts + strip together reproduce her block exactly,
                    for the sensitivity comparison.

Usage
  python build_count_areas.py            # writes data/count_areas_v6.gpkg
"""

import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd

HERE = Path(__file__).parent
POLYGONS = HERE / "data/source/eleanor_polygons_2026/PWRGC2026_CountBoundaries_Polygons2.shp"
LINES = HERE / "data/source/PW_RGC_2026_lines/PW RGC 2026.shp"
OUT = HERE / "data/count_areas_v6.gpkg"

BUFFER_M = 150        # Eleanor's buffer around transect counts
SNAP_M = 10           # tolerance for lines drawn on a polygon edge (1 km blocks)
UNNAMED = ("Danby", "Glaisdale")

# GWCT regions, as used in HEATHER_BEETLE_DAMAGE_SUMMER_2025.xlsx
GWCT_REGION = {
    **{m: "N.DALES" for m in ["Bollihope", "Crossgill", "Edmundbyers", "Eggleston", "Geltsdale",
                              "Grinton", "Moor House", "Stags Fell", "Wemmergill"]},
    **{m: "S.DALES/PEAK" for m in ["Dallowgill", "Ramsgill", "Swinton"]},
    **{m: "NYM" for m in ["Danby", "Rosedale", "Snilesworth", "Westerdale"]},
}


def log(msg):
    print(msg, flush=True)


def main():
    polys = gpd.read_file(POLYGONS).to_crs(27700)
    lines = gpd.read_file(LINES).to_crs(27700)
    lines["Moor"] = lines.Moor.fillna(UNNAMED[0])
    lines["Count"] = lines.Count.fillna(UNNAMED[1])

    # Assign each transect line to the polygon holding most of its length
    def home(line):
        inside = polys.geometry.buffer(SNAP_M).intersection(line).length / line.length
        return inside.idxmax() if inside.max() > 0.9 else None
    lines["fid"] = lines.geometry.apply(home)
    if lines.fid.isna().any():
        sys.exit(f"Lines outside every polygon:\n{lines[lines.fid.isna()][['Moor', 'Count']]}")

    names = lines.groupby("fid").apply(lambda g: sorted(set(zip(g.Moor, g.Count))))
    if missing := set(polys.index) - set(names.index):
        sys.exit(f"Polygons with no transects inside: {sorted(missing)}")

    rows = []
    for fid, keys in names.items():
        poly = polys.geometry[fid]
        if len(keys) == 1:
            moor, count = keys[0]
            rows.append(dict(Moor=moor, Count=count, Count_Type="Core", Boundary="Eleanor 2026",
                             src_fid=fid, geometry=poly))
            continue
        # One polygon, several count areas: rebuild each with the 150 m rule
        parts = []
        for moor, count in keys:
            sel = lines[(lines.Moor == moor) & (lines.Count == count)]
            g = sel.geometry.union_all().convex_hull.buffer(BUFFER_M).intersection(poly)
            parts.append(g)
            rows.append(dict(Moor=moor, Count=count, Count_Type="Core",
                             Boundary=f"Split from Eleanor 2026: transects + {BUFFER_M} m",
                             src_fid=fid, geometry=g))
        strip = poly.difference(gpd.GeoSeries(parts).union_all())
        moor = keys[0][0]
        rows.append(dict(Moor=moor, Count="Between counts", Count_Type="Strip",
                         Boundary="Eleanor 2026 minus the count areas above", src_fid=fid, geometry=strip))
        log(f"Split polygon {fid} ({poly.area / 1e4:,.0f} ha) into {len(keys)} {moor} count areas "
            f"+ {strip.area / 1e4:,.0f} ha of strip")

    out = gpd.GeoDataFrame(rows, crs=27700)
    out["GWCT_Region"] = out.Moor.map(GWCT_REGION)
    out["area_ha"] = (out.geometry.area / 1e4).round(2)
    out = out.sort_values(["Moor", "Count_Type", "Count"]).reset_index(drop=True)
    out = out[["Moor", "Count", "Count_Type", "GWCT_Region", "area_ha", "Boundary", "src_fid", "geometry"]]

    # Checks
    core = out[out.Count_Type == "Core"]
    expected = set(zip(lines.Moor, lines.Count))
    assert set(zip(core.Moor, core.Count)) == expected, "count areas do not match the named transects"
    assert not core.duplicated(["Moor", "Count"]).any(), "duplicate count areas"
    assert out.GWCT_Region.notna().all(), f"no GWCT region for {set(out[out.GWCT_Region.isna()].Moor)}"
    assert abs(out.area_ha.sum() - polys.area.sum() / 1e4) < 0.5, "total area differs from Eleanor's file"
    overlap = core.geometry.union_all().area
    assert abs(overlap - core.geometry.area.sum()) < 1, "count areas overlap"

    # How well does the 150 m rule reproduce Eleanor's own transect polygons?
    log("\nEleanor's polygon vs transects + 150 m, for the other transect counts (overlap / union):")
    for r in out[out.Boundary == "Eleanor 2026"].itertuples():
        sel = lines[(lines.Moor == r.Moor) & (lines.Count == r.Count)]
        if sel.geometry.union_all().distance(r.geometry.exterior) < SNAP_M:
            continue                                   # drawn block, not a buffered transect count
        rule = sel.geometry.union_all().convex_hull.buffer(BUFFER_M)
        iou = rule.intersection(r.geometry).area / rule.union(r.geometry).area
        log(f"  {r.Moor}, {r.Count}: {iou:.2f}")

    OUT.unlink(missing_ok=True)
    out.to_file(OUT, driver="GPKG", layer="count_areas")
    log(f"\nWrote {OUT}: {len(core)} count areas ({core.area_ha.sum():,.0f} ha), "
        f"{(out.Count_Type == 'Strip').sum()} strip row ({out[out.Count_Type == 'Strip'].area_ha.sum():,.0f} ha), "
        f"{core.Moor.nunique()} moors")
    with pd.option_context("display.width", 200, "display.max_rows", 100):
        log(out.drop(columns="geometry").to_string())


if __name__ == "__main__":
    main()
