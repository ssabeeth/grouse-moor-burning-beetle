---
title: "Satellite detection of heather burning and cutting on red grouse count areas"
subtitle: "Seasons 2017/18 to 2025/26. Version 6, October 2026. Draft for internal review, not for external circulation."
---

# Summary

We mapped ground that was **burnt or cut** each season (October to 15 April) on **28 red grouse count areas** (4,273 ha, 16 moors) in the North Pennines, Yorkshire Dales and North York Moors, for **nine seasons, 2017/18 to 2025/26**, using Sentinel-2 satellite imagery.

Burning and cutting are reported together. Both happen in the same season and cannot be told apart at 10 m resolution. Since the 2021 deep-peat regulation, estates on blanket bog have largely switched from burning to cutting, so the combined figure (**total heather management**) is the more useful quantity.

**Main findings**

- **On average 171 ha, 4.0% of the count area, was burnt or cut each season.** The lowest season was 2020/21 (2.5%), the highest 2024/25 (6.6%).
- **Management is sporadic at count-area level:** heavy in one season, little or nothing the next. Edmundbyers, Waskerley ran 37%, 0%, 12%, 0%, 41%, 0%, 1%, 40% and 2% across the nine seasons. This is how management now operates, so every season is reported.
- **The North York Moors are managed more intensively than the North Pennines:** 6.4% a year against 3.7%. This is the same ordering as Shewring et al. (2024), who report about 6.3% and 2.7% for burning alone.
- **Total management did not fall after the 2021 regulation:** 158 ha a season before it, 155 ha a season from 2021/22 to 2023/24. This fits the field observation that cutting replaced burning on blanket bog. It is compatible with Shewring's 73% fall in *burning*, which does not count cutting.
- **2024/25 is the highest season, 280 ha (6.6%).** It holds when the new Sentinel-2C satellite is left out (252 ha, still the highest) and at every threshold tested. It coincides with the heather beetle outbreak, consistent with PW's observation that managers burnt and cut beetle-killed heather (Section 3.5).
- **Negative controls behave as expected:** Moor House National Nature Reserve 0.1% a year, Geltsdale RSPB reserve 1.0%.
- **The findings are robust.** Changing the detection threshold by 0.02 either way shifts totals by 20 to 25% but leaves the order of seasons and moors, the regional contrast and the flat post-2021 total unchanged.
- **Small cuts are below what 10 m imagery can see.** Crossgill and Wemmergill show 0.2 to 0.3% a year, as PW anticipated: their cuts are mostly smaller than the 0.1 ha the method can resolve. Aerial photography would capture them.
- **New: Sentinel-2 can also estimate heather beetle damage.** Combining a simple rule with gradient boosting, trained without hand-drawn labels and tested on moors it had not seen, the estimate is within 0.2 of the July 2025 survey on 19 of 22 count areas and within 0.1 on 15 (Section 4).

**What changed since version 5**

This version follows an expert review of v5 by Phil (PW). His comments changed the method and, more importantly, how the results are read.

| PW's review | v5 | v6 |
|---|---|---|
| The season ends on 15 April; early April is often the busiest period | After-season images from 1 April | After-season images from **16 April** |
| Many fires and cuts are small, sometimes 10 m by 20 m | Minimum patch 0.5 ha | Minimum patch **0.1 ha** |
| Burning and cutting both happen October to 15 April | Tried to split them by timing | **One class:** burnt or cut |
| Every season is usable; management is sporadic | "Individual seasons must not be quoted" | Every season reported |
| 2024/25 reflects management of beetle-killed heather | 2024/25 excluded as anomalous | Included and examined |
| After 2021, blanket bog estates cut instead of burning | Expected a fall like Shewring's | Total management, so no fall is expected |
| Definitive boundaries are available | Boundaries reconstructed from lines | Eleanor's definitive polygons |
| Figure 1 did not show the proportion | Shaded polygons | Proportional circles and a count area by season table |

# 1. Background and objective

Prescribed burning on moorland managed for red grouse is contentious for its effects on peat, hydrology and habitat condition. Since 1 May 2021, burning on deep peat in protected sites in England has been banned unless licensed, and many estates on blanket bog now cut heather instead.

The objective is to measure, for each count area and season, **how much ground was burnt or cut**. Grouse are counted on these areas, so the results can later be related to the counts.

# 2. Methods

## 2.1 Count areas

The count area boundaries are the definitive polygons digitised by Eleanor (GWCT) in October 2026. Long-term counts are 1 km² blocks walked as six parallel transects 150 m apart. Other counts are two parallel transects, typically 3 km long and 500 m apart, with the count area drawn as the transects plus a 150 m buffer.

The polygon file carried no names, so each polygon was named from the count transects it contains (the 2026 transect file). This gives **28 count areas, 4,273 ha, on 16 moors**. One count area is new since v5: **Danby, Glaisdale**.

**Eggleston** was supplied as a single 1,928 ha block covering all seven of its count areas. We split it with the same rule (two transects plus 150 m), which reproduces Eleanor's other 13 transect polygons almost exactly (overlap 0.99 to 1.00). The 370 ha of strip between the Eggleston count areas is measured but left out of all totals. The choice makes no practical difference: Eggleston is 5.2% a year as seven count areas and 5.3% as the whole block; the North Pennines 3.7% against 3.9%.

Many count areas are considerably larger than the v5 reconstructions, so comparisons with v5 are made in percentages.

## 2.2 Detection

**Imagery.** Sentinel-2 Level-2A surface reflectance (10 m), with cloud, cloud shadow and snow masked using the s2cloudless cloud probability and the scene classification layer.

**Index.** The Normalised Burn Ratio, NBR = (NIR − SWIR2) / (NIR + SWIR2), from bands B8 and B12. It falls sharply when vegetation is burnt or cut.

**Two images per season**, each the median of all clear scenes in its window:

| Image | Window | Purpose |
|---|---|---|
| Before | 1 July to 30 September, before the season | Vegetation before any management that season |
| After | **16 April** to 30 June, after the season closes on 15 April | Ground burnt or cut during the season, including early April |

The drop in NBR between the two (dNBR) is the signal.

**Seasonal correction.** Heather browns over winter whether or not it is managed. Only a few per cent of the count areas are managed in any season, so the median dNBR across all count areas measures this background change. It is subtracted each season; afterwards the typical pixel shows no change (median between −0.004 and 0.004 in every season).

**Threshold.** A pixel counts as managed if its corrected dNBR exceeds **0.17**, the v5 calibration (the 95th percentile of dNBR in 2018/19). Re-checked on the new boundaries, the 2018/19 value is 0.176 against v5's 0.168, so the calibration holds. The 95th percentile varies between seasons (0.135 to 0.207) because the amount of management varies; matching the threshold to each season would force the same total every season, so one threshold is used throughout.

**Vegetation and land cover.** A pixel must have been vegetated in the before image (NDVI above 0.3), and trees, crops, buildings and water are masked using ESA WorldCover. Inside the count areas WorldCover classes 98.5% of the ground as grassland and none as shrubland, so this mask removes only the 1.5% of trees; it does not distinguish heather from bog.

**Minimum patch size: 0.1 ha** (ten connected pixels), the practical floor for 10 m imagery. Many fires and cuts are smaller than v5's 0.5 ha, and the smallest cuts (around 10 m by 20 m) are still below this floor.

**One class.** Burning and cutting are both carried out from October to 15 April, and what distinguishes a cut (straight edges, vehicle tracks a few metres wide) is not visible at 10 m. Results are therefore "burnt or cut".

**Image supply.** The after window had a median of 7.5 clear images per pixel. The lowest was Grinton, Long Gill in 2024/25, with about one clear image, so that value is less certain.

## 2.3 Robustness tests

- **Threshold:** the full analysis was repeated at 0.13, 0.15 and 0.19.
- **Sentinel-2C:** a third satellite entered service during 2024/25, raising the number of after-season scenes from about 106 to 156. 2024/25 was re-run without its scenes.
- **Eggleston boundaries:** seven count areas against Eleanor's whole block (Section 2.1).

# 3. Results

## 3.1 By season

![Ground burnt or cut within the 28 count areas, by season. The dashed line marks the deep-peat burning regulation of May 2021.](../outputs/figures/fig_season_totals.png){width=95%}

| Season | Burnt or cut (ha) | % of count area |
|---|---|---|
| 2017/18 | 169 | 4.0 |
| 2018/19 | 201 | 4.7 |
| 2019/20 | 156 | 3.7 |
| 2020/21 | 106 | 2.5 |
| 2021/22 | 199 | 4.7 |
| 2022/23 | 153 | 3.6 |
| 2023/24 | 114 | 2.7 |
| 2024/25 | 280 | 6.6 |
| 2025/26 | 158 | 3.7 |
| **Mean** | **171** | **4.0** |

## 3.2 By count area and moor

![Mean annual % of count area burnt or cut, by moor. Circle size and colour show the percentage.](../outputs/figures/fig_map.png){width=88%}

![Percentage of each count area burnt or cut, by season. Values of 0.5% or more are printed.](../outputs/figures/fig_heatmap.png){width=100%}

The count area by season figure shows the sporadic pattern clearly. Edmundbyers, Waskerley was heavily managed in four seasons (37%, 12%, 41% and 40%) and barely touched in the other five. Snilesworth, Lodge shows the same pattern. Most other working moors manage a few per cent of each count area in most seasons.

The imagery shows what a heavy season looks like on the ground. At Waskerley the managed ground in each heavy season is many small patches spread across the whole count area: the mosaic of rotational burning and cutting, not one event.

![Edmundbyers, Waskerley, season by season. Sentinel-2 false colour after each season; red is ground detected as burnt or cut; yellow is the count area.](../outputs_images/figures/contact_sheets/Edmundbyers_Waskerley.png){width=100%}

One single event stands out: **Bollihope, Pikestone Fell, 41% in 2022/23**, against 0% to 8% in every other season (v5 also found it). Unlike Waskerley's mosaic, it is one solid block of about 40 ha, far larger than a normal rotational burn of 2 to 5 ha. It may have been a wildfire, and is worth checking with the estate.

![Bollihope, Pikestone Fell, season by season. The 2022/23 detection is a single block of about 40 ha.](../outputs_images/figures/contact_sheets/Bollihope_PikestoneFell.png){width=100%}

![Mean annual % of count area burnt or cut, by moor and region.](../outputs/figures/fig_by_moor.png){width=80%}

Edmundbyers (14.7% a year) and Snilesworth (11.8%) are the most intensively managed moors. Five moors are at or below 1%: Geltsdale, Stags Fell, Crossgill, Wemmergill and Moor House.

## 3.3 Regional pattern

| Region | This study (burnt or cut) | Shewring et al. 2024 (burnt only) |
|---|---|---|
| North York Moors | 6.4% a year | about 6.3% a year |
| North Pennines | 3.7% a year | about 2.7% a year |
| Yorkshire Dales and Nidderdale | 3.7% a year | not compared |

The ordering matches Shewring et al., and the North York Moors figure is almost identical. The North Pennines figure is higher than theirs. Our figure includes cutting and theirs does not; with cutting replacing burning on blanket bog since 2021, a higher figure is to be expected. The two studies also cover different areas: count areas here, whole regions there.

\needspace{10\baselineskip}

## 3.4 The 2021 regulation

| Period | Mean burnt or cut per season |
|---|---|
| 2017/18 to 2020/21 (before) | 158 ha (3.7%) |
| 2021/22 to 2023/24 (after) | 155 ha (3.6%) |
| 2024/25 | 280 ha (6.6%) |
| 2025/26 | 158 ha (3.7%) |

Total management did not fall after the regulation. This is consistent with PW's field observation that estates on blanket bog switched from burning to cutting, which is not weather-limited, and with Shewring et al.'s 73% fall in burning, which excludes cutting. The satellite measures burning and cutting together, so it shows the total holding steady; the field observations show the switch within it.

## 3.5 2024/25 and the heather beetle outbreak

**2024/25 is the highest season in the record.** Without Sentinel-2C scenes it is 252 ha, still well above the next highest season (201 ha), and the ranking of count areas barely changes (rank correlation 0.97). It is the highest season at every threshold tested.

The peak coincides with the heather beetle outbreak and fits PW's field observation that managers burnt and cut beetle-killed heather to promote recovery. Count areas with more damage in the July 2025 survey were more heavily managed in 2024/25 (rank correlation 0.55). They were also the more heavily managed ones before the outbreak (0.70 in 2023/24), and relative to each count area's usual level the 2024/25 increase is spread across count areas, so the July 2025 survey on its own cannot single out the beetle response. Damage records from 2024, the summer before managers responded, would allow that comparison count area by count area.

Beetle-related management may also be under-counted here. Burning or cutting heather that is already dead changes the satellite signal less than burning live heather. This may be why the count areas with the heaviest damage show less management than usual in 2025/26 (rank correlation −0.40 relative to their usual level).

## 3.6 Sensitivity to the threshold

![Season totals at four thresholds. Levels shift; the pattern does not.](../outputs/report/fig_threshold_sensitivity.png){width=90%}

| Threshold | Total vs 0.17 | Season order* | Moor order* | NYM / NP ratio | After / before 2021 | Moor House | Geltsdale |
|---|---|---|---|---|---|---|---|
| 0.13 | +61% | 0.97 | 0.98 | 1.64 | 1.03 | 0.65% | 4.5% |
| 0.15 | +25% | 0.98 | 0.99 | 1.72 | 1.00 | 0.23% | 2.0% |
| **0.17** | main | 1 | 1 | 1.75 | 0.98 | 0.08% | 1.0% |
| 0.19 | −20% | 0.93 | 1.00 | 1.78 | 0.97 | 0.02% | 0.55% |

*\* Rank correlation with the 0.17 results.*

Between 0.15 and 0.19 the totals change by about a fifth, but every pattern in this report holds. At 0.13 the reserves start to register (Geltsdale 4.5% a year), so 0.13 is too permissive. Lowering the threshold barely changes the after / before 2021 ratio (0.97 to 1.03), so there is no sign that faint cut scars cluster after the regulation.

\needspace{14\baselineskip}

## 3.7 Checks

| Check | Expected | Result | Verdict |
|---|---|---|---|
| Moor House (National Nature Reserve) | Near zero | 0.08% a year | Pass |
| Geltsdale (RSPB reserve) | Low | 1.0% a year | Pass |
| North York Moors vs North Pennines | NYM higher (Shewring et al.) | 6.4% vs 3.7% | Pass |
| Crossgill and Wemmergill | More than Geltsdale (PW) | 0.26% and 0.23% vs 0.98% | Below detection: small cuts, as PW anticipated |
| 16 April window absorbs v5's "cut" class | v6 about equal to v5 burn + cut | 291 ha vs 211 + 65 = 276 ha (seven unchanged boundaries) | Pass |
| 2024/25 not caused by Sentinel-2C | Similar without it | 252 vs 280 ha, still the highest | Pass |

## 3.8 Comparison with version 5

On the 27 count areas common to both versions, the mean annual figure is almost unchanged: 4.10% in v5 and 4.03% in v6. Two things moved:

- **v5's "cut" class is now in the main result, as PW predicted.** On the seven count areas whose boundaries did not change, v5 found 211 ha of burning and 65 ha of later "cutting"; v6 finds 291 ha. This is consistent with PW's explanation: most of v5's "cutting" was early-April burning or cutting, diluted in v5's after image because that window began on 1 April, before the season closed.
- **Fewer exact zeros.** With the smaller minimum patch, 24% of count-area seasons show no management, against 39% in v5. The sporadic pattern remains, but small patches now register.

# 4. Heather beetle damage

## 4.1 Question

A heather beetle survey in July 2025 estimated the proportion of each count area damaged (23 count areas, 0 to 0.8). Can Sentinel-2 see the damage, and how close can it get to the survey's numbers?

## 4.2 Method

Beetle-killed heather turns red-brown in the summer it is attacked, then grey as the dead leaves drop. We built the estimate in two steps.

**Step 1: a rule-based estimate.** Two measures capture the damage:

- **Browning:** the drop in NBR from late summer 2024 to late summer 2025.
- **Greying:** the drop in colour saturation in the visible bands from early to late summer 2025 (flatter colour is greyer).

Summer 2025 was dry and the whole area browned, so the moor-wide median change is subtracted, as in the main method. A pixel is **flagged** if it browned by more than 0.10 or greyed by more than 0.05 beyond the moor-wide change. Ground burnt or cut in 2024/25 or 2025/26, plus a 20 m margin, is left out. The estimate is the **share of each count area flagged**. Nothing is fitted to the survey.

**Step 2: gradient boosting labelled automatically from the survey.** About 280 pixels were sampled in each count area, each described by its change in all ten Sentinel-2 bands (including the red-edge bands) and five indices across three time comparisons. Pixels from count areas the survey scored 0 were labelled healthy, and pixels from count areas scored 0.5 or more were labelled damaged; no labels were drawn by hand. A gradient boosting classifier learns from these labels, and each count area's mean probability is converted to a proportion. The model is always tested on a moor it has not seen (Section 4.4). The final estimate is the **average of the rule-based estimate and gradient boosting**. The formulae are in the Appendix.

## 4.3 Results

![Left: the rule-based estimate. Centre: the rule-based estimate averaged with gradient boosting. Right: the measures gradient boosting relies on most. Each moor is predicted without its own data; the grey band is within 0.2 of the survey.](../outputs/beetle_check/fig_beetle_gb.png){width=100%}

| | First attempt (raw browning) | Rule-based estimate | Rule-based + gradient boosting |
|---|---|---|---|
| Typical difference from the survey | 0.27 | 0.12 | **0.10** |
| Within 0.2 of the survey | | 19 of 22 | **19 of 22** |
| Within 0.1 of the survey | | 12 of 22 | **15 of 22** |
| Rank agreement with the survey | 0.58 | 0.81 | **0.88** |

*22 count areas: Grinton, Long Gill is left out as it had too few clear images. For comparison, predicting each count area from the average damage in its region, with no satellite data, gives a typical difference of 0.17.*

The two steps are complementary:

- **The rule captures the full range of damage.** Close matches include Danby, Low Moor (0.51 against 0.5), Eggleston, Ever Rigg (0.31 against 0.3) and Ramsgill, Raygill (0.22 against 0.2). Greying matters: Dallowgill, Bishops (survey 0.6) barely browned but greyed strongly.
- **Gradient boosting learns what is not damage.** At Moor House, Trout Beck the rule reads 0.55 against a survey score of 0, because the whole bog changed colour. Gradient boosting learned from the clean count areas that this greying is not damage and reads 0.05; the average reads 0.30.
- **Red-edge bands carry much of the signal.** After the spring-to-late-summer change in NBR, the measures gradient boosting relies on most are the red-edge bands B5 and B6, which respond to the loss of chlorophyll; four of its ten most-used measures are red-edge. Any future beetle model should include them.

**Other models.** We compared six models in Step 2: classic gradient boosting, gradient boosting in the LightGBM style, a Random Forest, Extra Trees, a support vector machine and logistic regression. Averaged with the rule, all gave typical differences between 0.10 and 0.12. Classic gradient boosting was best (0.098) and is used here. The consistent gain comes from combining a learned model with the rule; more and better labels will matter more than the choice of algorithm.

The damage is not clearly visible a year later: comparing summer 2026 with summer 2024 shows no clear relationship with the survey, because the heather recovered or because damaged heather was burnt or cut in 2025/26 (and is then masked out).

## 4.4 Checks for data leakage

Leakage means information about the answer reaching the model, which makes results look better than they will be on new data. The main burning and cutting analysis trains nothing on outcomes, so this applies to the beetle model. The design rules it out in these ways:

| Possible leak | How it is handled |
|---|---|
| A moor's own pixels or survey score used to predict it | Never: every estimate comes from a model trained without that moor (leave-one-moor-out) |
| Calibration fitted on the model's own training fit | No: the line converting probability to proportion uses out-of-fold probabilities, split by moor |
| Measures built from survey scores | No: all measures are image changes |
| Settings tuned on the test results | No: model settings were fixed in advance |
| Moor-wide reference includes the held-out moor | Imagery only, no survey scores; one moor barely moves a median over 4,000 ha |

Two further tests check what the design cannot:

- **Choosing the best model after the fact.** Picking the best of six models on the same 22 count areas flatters the winner. In a nested test the model was chosen inside each fold using the training moors only, then scored on the held-out moor. The combined estimate scored **0.105**, against 0.098, so the choice added very little; classic gradient boosting was chosen in 9 of the 15 folds.
- **Similar neighbouring moors.** When a whole region is hidden from training (the North York Moors, then the southern Dales), gradient boosting on its own does markedly worse on those 8 count areas (0.22 against 0.12): it benefits from having seen similar moors nearby. The combined estimate holds up (0.07, close to the rule's 0.08), because the rule needs no training. This is the main reason to report the combined estimate rather than the model alone.

## 4.5 Imagery where satellite and survey differ most

![The five count areas where the rule-based estimate and the survey differ most. Left to right: true colour 2024, 2025 and 2026; browning; greying; pixels counted as damage (red). Black is burnt or cut ground plus a 20 m margin.](../outputs/beetle_check/fig_evidence_disagree.png){width=100%}

At 10 m the true-colour images show the damage only faintly; the change maps say more. At these count areas the imagery mostly supports the survey:

- **Moor House, Trout Beck** (survey 0, rule 0.55, combined 0.30): the greying covers the whole bog, inside and outside the count area. It does not look like beetle damage, and gradient boosting discounts it.
- **Edmundbyers, Waskerley** (survey 0, rule 0.23): heavily burnt in 2024/25; the remaining flags sit around the burns.
- **Wemmergill, South Side** (survey 0.6, rule 0.37, combined 0.24): some flagged patches inside the count area, but fewer than the survey suggests; field photos would settle it.

![The five count areas where the rule-based estimate and the survey agree best, for comparison.](../outputs/beetle_check/fig_evidence_agree.png){width=100%}

## 4.6 Limits

- **Small sample, one year.** The flagging cut-offs and model settings were fixed in advance, but greying as a measure and the decision to average were settled while looking at the same 22 count areas; the nested test shows the model choice added little (Section 4.4). With 22 count areas, differences of 0.01 to 0.02 between methods are within chance. A second survey year would test the method on new data.
- **Clean ground reads slightly high.** Where the survey found little or no damage (0 to 0.1), the combined estimate runs about 0.11 high on average.
- **Survey scores are estimates in tenths,** so agreement much better than 0.1 cannot be expected.
- **Heather cannot yet be mapped separately from bog.** WorldCover has no heather class here, and an attempt to identify heather from its seasonal greenness did not help; learning from the clean count areas does part of this job.

# 5. Limits of the main analysis, and how to address them

- **Small cuts.** Cuts of a few hundred square metres are below the 0.1 ha floor, which is why Crossgill and Wemmergill show little. Aerial imagery would capture them.
- **Accuracy figure.** The threshold is calibrated on the data and the patterns are stable across every test, but an accuracy figure needs ground truth (Recommendation 1).
- **Burning and cutting together.** They cannot be separated at 10 m; aerial imagery would allow it.
- **Beetle damage and the signal.** Management of dead heather is likely under-counted; occasionally, damage appearing between the before and after images could be counted as management.
- **Image supply varies.** A few count-area seasons rest on very few clear images, notably Grinton, Long Gill in 2024/25.

# 6. Recommendations

**Immediate**

1. **Ground-truth the method.** Identify 15 to 20 sites where burning or cutting is known (from estate records or field visits), digitise them, and run the detection against them. This gives an accuracy figure and is the most important next step.
2. **Beetle damage evidence:** field photos from the July 2025 survey (ideally with locations), any damage scores from 2024, and a note of how the scores were estimated.
3. **Check two single results** with the estates: the 40 ha block at Pikestone Fell in 2022/23 (possibly a wildfire), and Geltsdale's detections (1.0% a year).

**Next**

4. **Aerial imagery (APGB, 25 cm or 50 cm, colour infrared).** At this resolution cuts and burns can be told apart by their shape, and small cuts become visible. Check whether GWCT already has access.
5. **Develop the beetle damage model.** Gradient boosting labelled automatically from the survey already improves the estimate (Section 4). A second survey year would test it on new data, and a few hundred hand-checked patches (beetle-brown, beetle-grey, burnt, cut, healthy; quickest from Google Earth's historical imagery) would let it map damage pixel by pixel. Test any version by holding out one moor at a time; a single 12/12 split of count areas leaves too few to train on.
6. **PlanetScope (3 m, near-daily),** free under Planet's Education and Research programme, is the intermediate option between Sentinel-2 and aerial imagery.
7. **Replace the fixed threshold with a supervised classifier** once ground-truth data exist.

\needspace{16\baselineskip}

# 7. Data and code

| File | Contents |
|---|---|
| `burn_pipeline.py` | Detection pipeline (Earth Engine Python API) |
| `build_count_areas.py` | Names Eleanor's polygons from the transect file; splits Eggleston |
| `beetle_check.py`, `beetle_rf.py`, `beetle_model_compare.py`, `beetle_evidence.py` | Beetle damage: rule-based estimate, pixel sampling, gradient boosting and model comparison with leakage tests, evidence sheets |
| `report_numbers.py` | Every number in this report, from the outputs below |
| `data/count_areas_v6.gpkg` | 28 count areas plus the Eggleston strip |
| `outputs/results_by_season.csv` | Results: 28 count areas by 9 seasons, threshold 0.17 |
| `outputs_t013/`, `outputs_t015/`, `outputs_t019/`, `outputs_no_s2c/` | Sensitivity runs |
| `outputs_images/` | Season-by-season image sheets for six example count areas |
| `outputs/beetle_check/` | Beetle damage comparison |

\newpage

# Appendix: formulae

Notation: $\rho_k(x,t)$ is the surface reflectance of Sentinel-2 band $k$ at pixel $x$ on date $t$ (digital numbers divided by 10,000). $\mathbf{1}[\cdot]$ is 1 when the condition holds and 0 otherwise.

## A.1 Indices

$$\mathrm{NBR} = \frac{\rho_{B8}-\rho_{B12}}{\rho_{B8}+\rho_{B12}}, \qquad \mathrm{NDVI} = \frac{\rho_{B8}-\rho_{B4}}{\rho_{B8}+\rho_{B4}}, \qquad \mathrm{NDMI} = \frac{\rho_{B8}-\rho_{B11}}{\rho_{B8}+\rho_{B11}}$$

$$\text{Brightness} = \tfrac{1}{3}\left(\rho_{B2}+\rho_{B3}+\rho_{B4}\right), \qquad \text{Saturation } S = \frac{\max(\rho_{B2},\rho_{B3},\rho_{B4})-\min(\rho_{B2},\rho_{B3},\rho_{B4})}{\max(\rho_{B2},\rho_{B3},\rho_{B4})}$$

$S$ is 0 for a perfectly grey surface and rises with colour.

## A.2 Composites

A scene is clear at $x$ if its s2cloudless cloud probability is below 40% and its scene classification is not cloud shadow, cloud (medium or high), cirrus or snow (classes 3, 8, 9, 10, 11). For any index $I$ and window $W$:

$$I_W(x) = \operatorname*{median}_{t \in W,\ \text{clear at } x} I(x,t)$$

## A.3 Change and seasonal correction (main analysis)

For season $s$, with the before window $W_b$ (1 July to 30 September) and the after window $W_a$ (16 April to 30 June):

$$\mathrm{dNBR}_s(x) = \mathrm{NBR}_{W_b}(x) - \mathrm{NBR}_{W_a}(x)$$

$$c_s = \operatorname*{median}_{x \in M} \mathrm{dNBR}_s(x), \qquad \mathrm{dNBR}^{*}_s(x) = \mathrm{dNBR}_s(x) - c_s$$

where $M$ is the moorland in all count areas, sampled at 100 m. The threshold check reports $p_{95,s}$, the 95th percentile of $\mathrm{dNBR}^{*}_s$ over $M$.

## A.4 Detection

$$r_s(x) = \mathbf{1}\!\left[\mathrm{dNBR}^{*}_s(x) > \tau\right]\cdot \mathbf{1}\!\left[\mathrm{NDVI}_{W_b}(x) > 0.3\right]\cdot \mathbf{1}\!\left[x \in M\right], \qquad \tau = 0.17$$

$$m_s(x) = r_s(x)\cdot \mathbf{1}\!\left[\,|P_s(x)| \ge 10\,\right]$$

where $P_s(x)$ is the patch of pixels with $r_s = 1$ connected to $x$ through any of its eight neighbours. Ten 10 m pixels are 0.1 ha.

## A.5 Areas and percentages

For count area $a$ with area $A_a$, and pixel area $\omega(x)$ (weighted by the fraction of the pixel inside $a$):

$$H_{a,s} = \sum_{x \in a} m_s(x)\,\omega(x), \qquad \%_{a,s} = 100\,\frac{H_{a,s}}{A_a}$$

For a group $G$ of count areas (a moor, a region, or all) over seasons $s = 1,\dots,9$, the mean annual percentage is area-weighted:

$$\%_G = 100\,\frac{\sum_{a\in G}\sum_s H_{a,s}}{\sum_{a\in G}\sum_s A_a}$$

## A.6 Beetle damage: rule-based estimate

$$\Delta\mathrm{NBR}(x) = \mathrm{NBR}_{\text{Jul-Sep 2025}}(x) - \mathrm{NBR}_{\text{Jul-Sep 2024}}(x), \qquad \Delta S(x) = S_{\text{Jul-Sep 2025}}(x) - S_{\text{16 Apr-30 Jun 2025}}(x)$$

$U$ is the usable ground: moorland, excluding pixels within 20 m of ground detected as burnt or cut in 2024/25 or 2025/26. Each change is centred on its median over $U$ in all count areas (sampled at 20 m):

$$\Delta\mathrm{NBR}^{*}(x) = \Delta\mathrm{NBR}(x) - \operatorname*{median}_{U}\Delta\mathrm{NBR}, \qquad \Delta S^{*}(x) = \Delta S(x) - \operatorname*{median}_{U}\Delta S$$

$$f(x) = \mathbf{1}\!\left[\Delta\mathrm{NBR}^{*}(x) < -0.10 \;\text{ or }\; \Delta S^{*}(x) < -0.05\right], \qquad R_a = \frac{1}{|a\cap U|}\sum_{x\in a\cap U} f(x)$$

## A.7 Beetle damage: gradient boosting

**Labels from the survey.** With $d_a$ the surveyed damage of count area $a$, each sampled pixel $x$ in $a$ gets
$$y_x = \begin{cases}0 & d_a = 0\\ 1 & d_a \ge 0.5\\ \text{no label} & \text{otherwise}\end{cases}$$

**Measures.** $z_x$ holds 45 changes: 15 measures (bands B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12; NBR, NDVI, NDMI, $S$, brightness) for each of three comparisons (summer 2025 minus summer 2024; late minus early summer 2025; summer 2026 minus summer 2024).

**Model.** Gradient boosting builds an additive score from $M$ small regression trees $h_m$ (depth 3):
$$F_M(z) = F_0 + \nu\sum_{m=1}^{M} h_m(z), \qquad p(z) = \frac{1}{1+e^{-F_M(z)}}$$
Each tree is fitted, on a random 80% of the labelled pixels, to the negative gradient of the log-loss
$$L = -\sum_x \left[\,y_x \log p(z_x) + (1-y_x)\log\!\left(1-p(z_x)\right)\right]$$
with learning rate $\nu = 0.05$ and $M = 150$ trees.

**From pixels to a count area.** The mean probability $\bar p_a = \frac{1}{n_a}\sum_{x\in a} p(z_x)$ is mapped to a proportion with a straight line fitted by least squares on the training count areas, using out-of-fold probabilities:
$$\hat d_a = \min\!\left(1,\ \max\!\left(0,\ \alpha + \beta\,\bar p_a\right)\right)$$

**Combined estimate.** $\hat c_a = \tfrac{1}{2}\left(R_a + \hat d_a\right)$.

## A.8 Evaluation

**Leave-one-moor-out.** For each moor $m$, the model and its calibration line are fitted using only the other moors' pixels and survey scores; $\hat c_a$ is then computed for the count areas of $m$. The moor-wide medians in A.3 and A.6 use imagery from all count areas but no survey scores.

$$\text{Typical difference} = \frac{1}{n}\sum_{a=1}^{n}\left|\hat c_a - d_a\right|, \qquad \text{Within } k = \#\left\{a : \left|\hat c_a - d_a\right| \le k\right\}$$

**Rank agreement** is Spearman's $\rho$: the Pearson correlation of the ranks of $\hat c_a$ and $d_a$.

**Paired bootstrap.** For $b = 1,\dots,5000$, draw $n$ count areas with replacement and compute $\Delta_b = \overline{|e^{\text{model}}|} - \overline{|e^{\text{rule}}|}$; the 90% interval is the 5th to 95th percentile of $\Delta_b$.

**Region-only benchmark.** $\hat d_a$ is the mean of $d$ over the other moors in the same GWCT region.

**Nested model choice.** For each held-out moor $m$, the training moors are split into five groups; each candidate model's combined estimate is scored by holding out each group in turn, the model with the lowest typical difference is refitted on all training moors, and it alone predicts $m$.

**Region hold-out.** As leave-one-moor-out, but every moor in one GWCT region is held out together (North York Moors; southern Dales). The northern Dales region cannot be held out: it contains every count area the survey scored 0.

# Reference

Shewring, M.P., Wilkinson, N.I., Teuten, E.L., Buchanan, G.M., Thompson, P. & Douglas, D.J.T. (2024) Annual extent of prescribed burning on moorland in Great Britain and overlap with ecosystem services. *Remote Sensing in Ecology and Conservation*. doi:10.1002/rse2.389
