# Heather burning, cutting and beetle damage on red grouse count areas

Sentinel-2 analysis, in Google Earth Engine, of heather management (burning plus cutting) on 28 red grouse count areas in northern England, seasons 2017/18 to 2025/26, and of heather beetle damage in summer 2025 against the July 2025 survey.

## Reports

| File | What it is |
|---|---|
| `report/report_v6.pdf` | Heather burning and cutting, version 2 (October 2026) |
| `report/beetle_report.pdf` | Heather beetle damage, summer 2025 |

Internal file names use v5 for version 1 and v6 for version 2.

## Code

| File | What it does |
|---|---|
| `build_count_areas.py` | Names Eleanor's count area polygons from the transect file; splits Eggleston (`data/count_areas_v6.gpkg`) |
| `burn_pipeline.py` | Detection of burnt or cut ground per count area and season |
| `v5_to_v6_changes.py` | Rebuilds version 1's boundaries and measures the effect of each change since version 1 |
| `beetle_check.py` | Beetle damage: image changes and the rule-based estimate, against the survey |
| `beetle_rf.py` | Samples pixels and their change measures for the beetle models |
| `beetle_model_compare.py` | Six models, gradient boosting, leakage tests |
| `beetle_unscored.py` | Beetle estimates for count areas the survey did not score |
| `beetle_evidence.py` | Image sheets where satellite and survey agree and differ |
| `report_numbers.py` | Every number quoted in the reports |
| `report/make_figures.py`, `report/make_beetle_figures.py` | Report figures and tables from the outputs |
| `report/build.sh` | Builds both PDFs (xelatex via latexmk) |

## Run

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
earthengine authenticate --scopes=https://www.googleapis.com/auth/earthengine,https://www.googleapis.com/auth/cloud-platform
export GEE_PROJECT=<your Earth Engine project>

python burn_pipeline.py --polygons data/count_areas_v6.gpkg            # main run, threshold 0.17
python beetle_check.py && python beetle_rf.py && python beetle_model_compare.py --nested
python beetle_unscored.py && python report_numbers.py
python report/make_figures.py && python report/make_beetle_figures.py && report/build.sh
```

Sensitivity runs: `--threshold 0.13|0.15|0.19 --tag t013` etc., `--seasons 2024 --exclude-s2c --tag no_s2c`. Always pass `--tag` with `--only` so `outputs/` is not overwritten.
