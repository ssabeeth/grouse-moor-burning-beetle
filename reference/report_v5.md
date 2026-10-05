# Satellite Detection of Burning and Cutting on Red Grouse Count Areas

**Extent of moorland scarring, 2017/18 – 2024/25**

*Draft — internal review. Not for external circulation.*

---

## Summary

We mapped burnt and cut moorland across **27 red grouse count areas** (2,537 ha, 16 moors in the North Pennines, North York Moors and Yorkshire Dales) using Sentinel-2 satellite imagery in Google Earth Engine, across **eight burning seasons**.

**What we can stand behind:**

- The **spatial pattern is corroborated** against published work. Our North York Moors figure (8.6%/yr) and North Pennines figure (3.3%/yr) reproduce both the magnitude and the ~2.5× ratio reported by Shewring et al. (2024), derived independently by a different method.
- **Edmundbyers (Waskerley), Snilesworth (Lodge) and Rosedale (Pits)** are the most heavily scarred count areas, at 11–15% per year.
- **Moor House (0.03%/yr) and Geltsdale (0.4%/yr)** show almost no scarring — as expected for a National Nature Reserve and an RSPB reserve. This is a useful negative control.

**What we cannot stand behind, and why it matters:**

- **Individual season figures must not be quoted.** Per-season detection behaves in an all-or-nothing manner: a count area registers either a large fraction of its area or exactly zero, with nothing in between (§4.3). This is not how rotational burning behaves, and it indicates a threshold instability in our method.
- **The between-year trend is unreliable.** We do *not* reproduce the 73% fall in English burning after the May 2021 deep-peat regulation that Shewring et al. report. Our data show an apparent *rise*. This reflects a fault in our method, not a fact about the world.
- **2024/25 is anomalous** (2.5× any other season) and coincides with a step change in satellite image supply. It is excluded from interpretation.
- **Cutting cannot be separated from burning** at 10 m resolution. Report the combined figure.

> **In short: the eight-season mean per count area is usable. Any single season is not.**

---

## 1. Background and objective

Prescribed burning on moorland managed for red grouse is contentious for its effects on peat soils, hydrology and habitat condition. Legislation introduced in England on 1 May 2021 prohibits burning on deep peat within protected areas.

The objective was to quantify, per count area and per season: **how much ground was burnt or cut, and to what extent**.

An initial approach using Google Earth Pro was abandoned. At the zoom levels needed to cover 16 moors, burn scars (typically 30–100 m wide) fall below the resolvable scale, and the historical imagery is a cloud-free composite rather than date-specific coverage — so before/after comparison is not possible.

---

## 2. Methods

### 2.1 Count area boundaries — an important caveat

The supplied shapefile (`PW_RGC_2026`) contains **35 POLYLINE features, not polygons** — 132 km of lines with no enclosed area. Burn *extent* cannot be derived from lines.

Inspection showed the lines were drawn inconsistently, under three different conventions:

| Convention | Records |
|---|---|
| Closed ring (a polygon in all but name) | 3 |
| Near-closed outline (traced, never snapped shut — gaps of 4–9 m) | 4 |
| Pair or set of edge lines bounding a block | 20 |
| Parallel transect set (Stags Fell ×6, Geltsdale ×4) | 8 |

We reconstructed **27 count area polygons (2,537 ha)** using the closed ring where one existed (preserving concavity) and the convex hull of the line set otherwise.

> **These polygons are a reconstruction, not an authoritative boundary.** Eggleston alone accounts for 872 ha — over a third of the total — so its boundaries carry disproportionate weight in any comparison between moors. **If authoritative count area boundaries exist within the organisation, they should replace these.** This is the single largest source of uncertainty in the study and it is not a remote-sensing problem.

![Figure 1](fig1_map.png)

### 2.2 Detection

Following the broad approach of Shewring et al. (2024), with modifications:

**Imagery.** Sentinel-2 Level-2A surface reflectance (10 m), cloud- and snow-masked using the s2cloudless probability layer and the Scene Classification Layer.

**Index.** The Normalised Burn Ratio, which falls sharply when vegetation is removed:

> NBR = (NIR − SWIR₂) / (NIR + SWIR₂)  — Sentinel-2 bands B8 and B12

**Three composites per season:**

| Composite | Window | Purpose |
|---|---|---|
| Baseline | Jul–Sep, pre-season | Vegetation before any burning that year |
| Burn window | Apr–May, post-season | Scar present by the close of the legal burn window? |
| Late summer | Aug–Sep, post-season | Scar appeared only afterwards? |

**Seasonal correction.** Heather senesces in winter; NBR falls moor-wide whether or not anything burns. Because only a few percent of moorland burns in any season, the **median dNBR across all count areas is the seasonal background**. We subtract it per composite. Verification: the corrected median came out at **−0.004**, i.e. the typical moorland pixel shows no change — as it should.

**Threshold — measured, not assumed.** The dNBR distribution across the count areas for 2018/19:

| Percentile | p50 | p90 | p95 | p97 | p99 | max |
|---|---|---|---|---|---|---|
| dNBR | −0.004 | 0.113 | **0.168** | 0.223 | 0.340 | 0.473 |

The long right tail is the burning signal. We set the threshold at **0.17 (p95)**, implying ~5% of moorland disturbed in a heavy season — consistent with published figures for these moors.

**Moorland mask.** ESA WorldCover (10 m): shrubland, grassland, herbaceous wetland and moss/lichen retained; cropland, built-up, woodland and water excluded. **This is essential** — silage cutting and harvest on enclosed farmland are spectrally near-identical to burning, and without the mask the classifier fires across the entire agricultural landscape.

**Minimum patch size.** 0.5 ha. Recommended burn sizes are 2–5 ha (DEFRA), so smaller detections are noise. Applied to each class independently (see §4.5).

![Figure 6](fig6_detection_eggleston.png)
*Detected scarring (red) within an Eggleston count area, 2018/19. Blue = scar appearing only after the burn window closed — largely an artefact, see §4.5.*

![Figure 7](fig7_falsecolour.png)
*The same area in false colour (SWIR–NIR–red), where burn scars are far more legible than in true colour.*

---

## 3. Results

### 3.1 Overall

**831 ha of scarring detected across eight seasons**, a mean of **104 ha/yr — 4.1% of the count area estate per year**.

### 3.2 Spatial pattern

![Figure 3](fig3_bymoor.png)

| Rank | Count area | Mean annual % scarred |
|---|---|---|
| 1 | Edmundbyers — Waskerley | 15.5% |
| 2 | Snilesworth — Lodge | 11.9% |
| 3 | Rosedale — Pits | 10.7% |
| 4 | Danby — Low Moor | 8.5% |
| 5 | Bollihope | 6.7% |
| … | | |
| 15 | Crossgill | 0.1% |
| 16 | Moor House | **0.03%** |

### 3.3 Corroboration against published work

![Figure 4](fig4_validation.png)

| Region | This study | Shewring et al. (2024) |
|---|---|---|
| North York Moors | **8.6%/yr** | ~6.3%/yr |
| North Pennines | **3.3%/yr** | ~2.7%/yr |

Same ordering, comparable magnitude, and near-identical ratio between regions — from an independent method. **Moor House at 0.03%** is a strong negative control: it is a National Nature Reserve, not a working grouse moor, and the method correctly finds nothing there.

---

## 4. Discussion — where this fails

### 4.1 The regulation effect is not reproduced

![Figure 2](fig2_trend.png)

| Period | Mean scarring |
|---|---|
| 2017/18 – 2020/21 (pre-regulation) | 3.15%/yr |
| 2021/22 – 2024/25 (post-regulation) | **5.05%/yr** |

Shewring et al. found English burning fell **73%** in 2021/22 following the deep-peat regulation. **We show the opposite.**

This was our principal validation test and **it failed.** Given that the published result rests on a validated deep neural network (84.9% balanced accuracy) applied nationally, and ours on a thresholded two-composite classifier applied to 27 reconstructed polygons, **the most likely explanation is a fault in our method.**

**No claim about the effect of the regulation should be made on this evidence.**

### 4.2 Detection is all-or-nothing between seasons — the clearest evidence of a fault

Examining a single count area across all eight seasons makes the problem unmistakable.

![Figure 8](fig8_waskerley_timeseries.png)

**Edmundbyers – Waskerley (66.8 ha), our most heavily burnt count area:**

| Season | Detected | % of count area | Usable images |
|---|---|---|---|
| 2017/18 | 19.3 ha | 28.9% | 32 |
| 2018/19 | **0.0 ha** | 0.0% | 33 |
| 2019/20 | 7.7 ha | 11.6% | 37 |
| 2020/21 | **0.0 ha** | 0.0% | 28 |
| 2021/22 | 22.8 ha | 34.3% | 35 |
| 2022/23 | **0.0 ha** | 0.0% | 25 |
| 2023/24 | **0.0 ha** | 0.0% | 26 |
| 2024/25 | 32.5 ha | **48.9%** | 41 |

The pattern is **binary**: either a third of the area, or exactly nothing. There are no intermediate values in eight years.

**Rotational burning does not behave like this.** Burning proceeds on 8–25 year rotations, so a working moor should show a few percent burnt in most years — a sequence like 2%, 5%, 1%, 4%. Alternating between a third of the moor and zero is not a management regime; it is an artefact.

Three things rule out the innocent explanations:

- **It is not cloud.** The zero seasons have 25–37 usable images per pixel — ample.
- **48.9% in one season is not credible.** Burning half a count area in a single year would be extraordinary and would be visible on the ground.
- **The zeros are exact.** Not 0.3 ha or 0.1 ha — zero. A genuine low-burning year would still register scattered small patches.

**Most likely mechanism.** The seasonal correction (§2.2) subtracts a single moor-wide median dNBR from every pixel, recomputed each season. If a count area's dNBR distribution sits close to the 0.17 threshold, a small year-to-year shift in that offset can push the *whole area* above or below the threshold at once. That would produce exactly this binary behaviour — and would also explain why our between-year trend contradicts the published one.

The correction is necessary (without it, winter senescence swamps the signal — see §2.2), but applying it as a **single global offset** appears to be too blunt. A per-count-area or per-scene offset, or a supervised classifier that does not depend on a hard threshold at all, would likely resolve this.

### 4.3 What this means for the results

The **eight-season mean per count area remains usable**, because the all-or-nothing behaviour averages out across seasons — and because the resulting spatial pattern independently reproduces the published figures (§3.3) and correctly identifies the non-grouse-moor sites as near-zero.

**Individual season figures should not be quoted, and no temporal claim should be made.**

### 4.4 2024/25 is anomalous and is excluded

![Figure 5](fig5_anomaly.png)

2024/25 shows **227.8 ha** scarred, against a 104 ha eight-season mean — 2.5× any other season, with 14 of 27 count areas up.

It also shows a sharp jump in image supply: **39 usable images per pixel, against 19–33 in every other season.**

**Sentinel-2C launched in September 2024** and entered service during this season. A third satellite yields more images, and potentially small inter-sensor calibration differences. A systematic shift entering the archive at precisely the moment our anomaly appears is unlikely to be coincidence, though we have not confirmed it.

**2024/25 is excluded from interpretation pending investigation.**

### 4.5 Cutting cannot be separated from burning

Our design attempted to separate the two by *timing*: burning is legally restricted to roughly 1 October – 15 April, whereas cutting is not. A scar present by the close of the burn window is likely a burn; one appearing only over summer is likely a cut.

**This does not work reliably.**

An initial run produced blue "cut" pixels forming a **thin fringe around every red burn patch** — an artefact of applying the noise filter to the combined scar and splitting afterwards, so single edge pixels survived. Fixing this (filtering each class independently) removed most of it, but a residual remains that we cannot attribute with confidence.

The fundamental limit is **resolution, not spectra**. What makes cutting recognisable is *geometry* — straight edges, sharp corners, vehicle tracks. Vehicle tracks are 2–3 m wide: **physically invisible** at 10 m. This is why Shewring et al. also could not separate them, finding 28.1% of known cuts classified as burns.

> **Report the combined "burnt or cut" figure.** Do not report the split.

---

## 5. Recommendations

### Immediate

1. **Obtain authoritative count area boundaries.** The largest uncertainty in this study is not the satellite method — it is that the boundaries are reconstructed from inconsistently-drawn lines. This is fixable within the organisation, today.
2. **Ground-truth the classifier.** Digitise 15–20 polygons in Google Earth Pro over sites where burning or cutting can be confidently identified by eye, and run the classifier against them. This produces a **confusion matrix and an accuracy figure** — the difference between "here are some numbers" and "here are some numbers, validated at X% accuracy."
3. **Fix the threshold instability (§4.2).** The single global seasonal offset is the prime suspect for the all-or-nothing behaviour. Test a per-count-area offset, and test whether the binary pattern persists. Until this is resolved, no season-level figure can be quoted.
4. **Investigate the Sentinel-2C hypothesis** by re-running 2024/25 with 2C scenes excluded.

### Next

5. **PlanetScope (3 m, 4-band, near-daily).** Free under Planet's Education & Research programme. At 3 m a cut strip is 10 pixels across rather than 3 — geometry becomes legible *while the spectral maths is retained*. This is the single largest available upgrade, and is what Shewring et al. themselves recommend.
6. **Replace the threshold with a supervised classifier** (Random Forest in Earth Engine). Trained on the ground-truth points from (2), this removes the hard threshold entirely — which would address §4.2 at root rather than patching it — and yields a quotable accuracy figure.
7. **APGB aerial imagery (25 cm, colour-infrared).** Check whether the organisation already holds a licence — many conservation bodies do. This is the resolution at which cutting-vehicle tracks appear, and the only route to a defensible cut/burn split.

---

## 6. Data and code

| File | Contents |
|---|---|
| `gee_burn_v5.js` | Detection script (Earth Engine, self-contained) |
| `gee_export_images.js` | Batch export of before/after/detection panels per count area |
| `gee_burn_animations.js` | 8-season animations per count area (export to Drive; the inline Console preview rate-limits and should not be relied on) |
| `gee_diagnostic.js` | Threshold percentiles, image counts, per-step pixel loss |
| `PW_RGC_2026_count_areas.geojson` | The 27 reconstructed count area polygons |
| `grouse_burn_v5.csv` | Results: 27 count areas × 8 seasons |

---

## Reference

Shewring, M.P., Wilkinson, N.I., Teuten, E.L., Buchanan, G.M., Thompson, P. & Douglas, D.J.T. (2024) Annual extent of prescribed burning on moorland in Great Britain and overlap with ecosystem services. *Remote Sensing in Ecology and Conservation*. doi:10.1002/rse2.389
