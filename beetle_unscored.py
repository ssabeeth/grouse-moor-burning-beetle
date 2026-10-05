#!/usr/bin/env python3
"""
Beetle damage estimates for count areas the July 2025 survey did not score (five at Eggleston).

Same method as the main estimate in beetle_model_compare.py: gradient boosting trained on the
automatically labelled pixels of all 22 scored count areas with enough images, its probability
mapped to a proportion with a straight line fitted on out-of-fold probabilities (split by moor),
then averaged with the rule-based estimate from beetle_check.py. Uses the pixels already cached
in rf_samples.csv, so no Earth Engine call.

Usage
  python beetle_unscored.py      # writes outputs/beetle_check/unscored_estimates.csv
"""
from pathlib import Path

import pandas as pd

from beetle_check import COMPARISONS
from beetle_model_compare import MAIN, cv_estimates, make_models
from burn_pipeline import STRIP, load_polygons, log

OUT = Path("outputs/beetle_check")


def main():
    core = load_polygons("data/count_areas_v6.gpkg")
    core = core[core.Count_Type != STRIP]
    check = pd.read_csv(OUT / "beetle_check.csv")
    survey = pd.read_csv(OUT / "beetle_vs_survey.csv")
    scored = survey[survey.enough_images][["Moor", "Count", "damage"]]
    todo = check[check.damage.isna() & check.enough_images][["Moor", "Count", "share_flagged__v1_edges"]]

    px = pd.read_csv(OUT / "rf_samples.csv").merge(core[["uid", "Moor", "Count", "GWCT_Region"]], on="uid")
    px = px.merge(scored, on=["Moor", "Count"], how="left")
    keys = set(zip(todo.Moor, todo.Count))
    px["group"] = ["unscored" if (m, c) in keys else "scored" for m, c in zip(px.Moor, px.Count)]
    px = px[(px.group == "unscored") | px.damage.notna()].reset_index(drop=True)
    cols = [c for c in px.columns if any(c.endswith(f"__{k}") for k in COMPARISONS)]

    gb = cv_estimates(px, todo, cols, make_models()[MAIN], [["unscored"]], by="group")
    res = todo.rename(columns={"share_flagged__v1_edges": "rule"}).assign(gradient_boosting=gb)
    res["combined"] = (res.rule + res.gradient_boosting) / 2
    res.round(3).to_csv(OUT / "unscored_estimates.csv", index=False)
    log(f"{len(res)} unscored count areas, trained on {scored.shape[0]} scored ones:\n"
        + res.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
