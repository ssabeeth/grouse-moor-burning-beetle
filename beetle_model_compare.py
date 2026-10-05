#!/usr/bin/env python3
"""
Beetle damage models: compares six models, makes gradient boosting the main one, and tests for
leakage. Uses the pixels cached by beetle_rf.py (run that first).

Design, for every model: labels taken automatically from the July 2025 survey (survey-0 count
areas healthy, >= 0.5 damaged; others unlabelled), change measures only, leave-one-moor-out. A
count area's estimate is its mean pixel probability, mapped to a proportion with a straight line
fitted on the training count areas, using out-of-fold probabilities (5 folds split by moor) so the
line is never fitted to a model's own training fit. Each model is also averaged with the
rule-based estimate from beetle_check.py.

Leakage tests
  bootstrap   paired bootstrap over count areas: 90% interval of each model's typical error minus
              the rule's. An interval spanning zero means the gap is within chance.
  region      the main model predicting a whole region it never saw (North York Moors; southern
              Dales). N.DALES cannot be held out: it holds every survey-0 count area.
  nested      (--nested, about 15 minutes on 10 cores) the "pick the best of six" procedure scored
              honestly: within each held-out moor's fold, the model is chosen using the training
              moors only (5-fold split by moor), then the held-out moor is predicted with it.

Usage
  python beetle_model_compare.py              # comparison, main model outputs, region test
  python beetle_model_compare.py --nested     # also the nested selection test
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import clone
from sklearn.ensemble import (ExtraTreesClassifier, GradientBoostingClassifier,
                              HistGradientBoostingClassifier, RandomForestClassifier)
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from beetle_check import COMPARISONS
from beetle_rf import REGION_COLOUR, area_means, score
from burn_pipeline import STRIP, load_polygons, log

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = Path("outputs/beetle_check")
MAIN = "Gradient boosting (classic)"


def make_models(n_jobs=-1):
    return {
        "Random Forest": RandomForestClassifier(n_estimators=300, min_samples_leaf=20, max_features="sqrt",
                                                n_jobs=n_jobs, random_state=0),
        "Extra Trees": ExtraTreesClassifier(n_estimators=300, min_samples_leaf=20, max_features="sqrt",
                                            n_jobs=n_jobs, random_state=0),
        "Gradient boosting (histogram, LightGBM-style)": HistGradientBoostingClassifier(
            learning_rate=0.05, max_iter=200, max_leaf_nodes=15, min_samples_leaf=40, l2_regularization=1.0,
            random_state=0),
        "Gradient boosting (classic)": GradientBoostingClassifier(n_estimators=150, learning_rate=0.05,
                                                                  max_depth=3, subsample=0.8, random_state=0),
        "SVM (RBF kernel)": make_pipeline(StandardScaler(), SVC(C=1.0, probability=True, random_state=0)),
        "Logistic regression": make_pipeline(StandardScaler(), LogisticRegression(C=0.1, max_iter=2000)),
    }


def cv_estimates(px, areas, cols, model, holdouts, by="Moor"):
    """Estimates for the count areas in each held-out group, from a model trained without that group.
    `holdouts` is a list of lists of group values (moors, or regions with by="GWCT_Region")."""
    est = {}
    for group in holdouts:
        tr, te = px[~px[by].isin(group)], px[px[by].isin(group)]
        lab = tr[(tr.damage == 0) | (tr.damage >= 0.5)]
        y = (lab.damage >= 0.5).values
        folds = GroupKFold(min(5, lab.Moor.nunique()))
        p_tr = pd.Series(np.nan, index=tr.index)
        p_tr.loc[lab.index] = cross_val_predict(clone(model), lab[cols], y, groups=lab.Moor, cv=folds,
                                                method="predict_proba")[:, 1]
        m = clone(model).fit(lab[cols], y)
        rest = tr.index.difference(lab.index)
        if len(rest):
            p_tr.loc[rest] = m.predict_proba(tr.loc[rest, cols])[:, 1]
        cal = area_means(tr, p_tr.values).to_frame("p").join(tr.groupby(["Moor", "Count"]).damage.first())
        slope, icpt = np.polyfit(cal.p, cal.damage, 1)
        est.update((icpt + slope * area_means(te, m.predict_proba(te[cols])[:, 1])).to_dict())
    return areas.set_index(["Moor", "Count"]).index.map(est).to_numpy(dtype=float).clip(0, 1)


def lomo(px, areas, cols, model):
    return cv_estimates(px, areas, cols, model, [[m] for m in areas.Moor.unique()])


def boot_diff(err_a, err_b, n=5000, seed=0):
    """90% interval of mean(err_a) - mean(err_b), resampling count areas."""
    idx = np.random.default_rng(seed).integers(0, len(err_a), size=(n, len(err_a)))
    return np.percentile(err_a[idx].mean(axis=1) - err_b[idx].mean(axis=1), [5, 95]).round(3)


def nested_fold(px, areas, cols, moor):
    """Choose the model using the training moors only, then predict the held-out moor with it."""
    models = make_models(n_jobs=1)
    tr_px, tr_areas = px[px.Moor != moor], areas[areas.Moor != moor].reset_index(drop=True)
    moors = np.random.default_rng(0).permutation(sorted(tr_areas.Moor.unique()))
    inner = [list(g) for g in np.array_split(moors, 5)]
    errs = {}
    for name, model in models.items():
        e = cv_estimates(tr_px, tr_areas, cols, model, inner)
        errs[name] = np.abs((e + tr_areas.satellite_estimate.values) / 2 - tr_areas.damage.values).mean()
    best = min(errs, key=errs.get)
    est = cv_estimates(px, areas, cols, models[best], [[moor]])
    return moor, best, est


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nested", action="store_true")
    a = ap.parse_args()

    core = load_polygons("data/count_areas_v6.gpkg")
    core = core[core.Count_Type != STRIP]
    survey = pd.read_csv(OUT / "beetle_vs_survey.csv")
    areas = survey[survey.enough_images].reset_index(drop=True)
    px = (pd.read_csv(OUT / "rf_samples.csv").merge(core[["uid", "Moor", "Count", "GWCT_Region"]], on="uid")
          .merge(areas[["Moor", "Count", "damage"]], on=["Moor", "Count"]).dropna().reset_index(drop=True))
    cols = [c for c in px.columns if any(c.endswith(f"__{k}") for k in COMPARISONS)]
    y, rule = areas.damage.values, areas.satellite_estimate.values
    err_rule = np.abs(rule - y)
    log(f"{len(px):,} pixels, {len(cols)} change measures, {len(areas)} count areas")

    # 1. Six models, leave-one-moor-out
    rows = [{"model": "Rule-based estimate (no model)", "combined_with_rule": False, **score(areas, rule),
             "vs_rule_90pc": ""}]
    preds = areas[["GWCT_Region", "Moor", "Count", "damage"]].assign(rule=rule)
    for name, model in make_models().items():
        log(f"  {name} ...")
        est = lomo(px, areas, cols, model)
        for combined, e in [(False, est), (True, (est + rule) / 2)]:
            lo, hi = boot_diff(np.abs(e - y), err_rule)
            rows.append({"model": name, "combined_with_rule": combined, **score(areas, e),
                         "vs_rule_90pc": f"{lo:+.3f} to {hi:+.3f}"})
        preds[name] = est.round(3)
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "model_comparison.csv", index=False)
    preds["main_combined"] = ((preds[MAIN] + preds.rule) / 2).round(3)
    preds.to_csv(OUT / "model_comparison_predictions.csv", index=False)
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        log("\n" + res.to_string(index=False))

    # 2. Main model: which measures it uses (fit on all labelled pixels)
    lab = px[(px.damage == 0) | (px.damage >= 0.5)]
    main = clone(make_models()[MAIN]).fit(lab[cols], lab.damage >= 0.5)
    imp = pd.Series(main.feature_importances_, index=cols).sort_values(ascending=False)
    imp.round(4).to_csv(OUT / "gb_importance.csv", header=["importance"])
    log("\nTop measures (" + MAIN + "):\n" + imp.head(10).round(3).to_string())

    # 3. Region hold-out
    reg = areas[["GWCT_Region", "Moor", "Count", "damage"]].assign(lomo=preds[MAIN].values, rule=rule)
    reg["region_out"] = np.nan
    for region in ["NYM", "S.DALES/PEAK"]:
        e = cv_estimates(px, areas, cols, make_models()[MAIN], [[region]], by="GWCT_Region")
        sel = (areas.GWCT_Region == region).values
        reg.loc[sel, "region_out"] = e[sel]
    held = reg.dropna(subset=["region_out"])
    region_rows = []
    for label, col in [("rule", "rule"), ("main model, leave-one-moor-out", "lomo"),
                       ("main model, whole region held out", "region_out")]:
        est = held[col].values if col != "lomo" else held.lomo.values
        comb = (est + held.rule.values) / 2 if col != "rule" else est
        region_rows.append({"estimate": label, "typical_error": round(np.abs(est - held.damage).mean(), 3),
                            "combined_typical_error": round(np.abs(comb - held.damage).mean(), 3)})
    rr = pd.DataFrame(region_rows)
    rr.to_csv(OUT / "region_holdout.csv", index=False)
    held.round(3).to_csv(OUT / "region_holdout_areas.csv", index=False)
    log(f"\nRegion hold-out ({len(held)} count areas in NYM and S.DALES/PEAK):\n" + rr.to_string(index=False))

    # 4. Figure: rule vs rule + main model, and importances
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), gridspec_kw={"width_ratios": [1, 1, 0.9]})
    for ax, (col, title) in zip(axes[:2], [("rule", "Rule-based estimate (nothing fitted)"),
                                           ("main_combined", "Rule-based estimate averaged with gradient boosting")]):
        s = score(areas, preds[col].values)
        ax.plot([0, 1], [0, 1], color="#999", lw=0.9, zorder=1)
        ax.fill_between([0, 1], [-0.2, 0.8], [0.2, 1.2], color="#eee", zorder=0)
        for region, g in preds.groupby("GWCT_Region"):
            ax.scatter(g.damage, g[col], s=40, color=REGION_COLOUR.get(region, "#888"), edgecolor="white",
                       lw=1.0, label=region, zorder=3)
        for _, row in preds.iterrows():
            ax.annotate(row.Count, (row.damage, row[col]), xytext=(4, 3), textcoords="offset points",
                        fontsize=6.5, color="#444")
        ax.set_title(f"{title}\ntypical error {s['typical_error']:.2f}, within 0.2: {s['within_0.2']}/{len(preds)}, "
                     f"within 0.1: {s['within_0.1']}/{len(preds)}, rank r {s['rank_r']:.2f}",
                     loc="left", fontsize=9.5, weight="bold")
        ax.set_xlabel("Surveyed beetle damage, July 2025")
        ax.set_xlim(-0.05, 1.0)
        ax.set_ylim(-0.05, 1.0)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Satellite estimate (each moor predicted without its own data)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    top = imp.head(12)[::-1]
    axes[2].barh([c.replace("__", ", ") for c in top.index], top.values, color="#2a78d6")
    axes[2].set_title("Measures gradient boosting relies on most", loc="left", fontsize=9.5, weight="bold")
    axes[2].tick_params(labelsize=8)
    axes[2].spines[["top", "right"]].set_visible(False)
    fig.suptitle("Beetle damage: gradient boosting labelled automatically from the survey, tested one moor at a "
                 "time. Grey band = within 0.2", x=0.01, ha="left", fontsize=12, weight="bold")
    plt.tight_layout()
    fig.savefig(OUT / "fig_beetle_gb.png", dpi=160, facecolor="white")

    # 5. Nested selection test
    if a.nested:
        log("\nNested selection test (model chosen inside each fold from the training moors only) ...")
        out = Parallel(n_jobs=-1, verbose=5)(delayed(nested_fold)(px, areas, cols, m) for m in areas.Moor.unique())
        est, chosen = np.full(len(areas), np.nan), {}
        for moor, best, e in out:
            sel = (areas.Moor == moor).values
            est[sel] = e[sel]
            chosen[moor] = best
        comb = (est + rule) / 2
        lo, hi = boot_diff(np.abs(comb - y), err_rule)
        nested = pd.DataFrame([{"procedure": "nested: best of six chosen inside each fold, averaged with rule",
                                **score(areas, comb), "vs_rule_90pc": f"{lo:+.3f} to {hi:+.3f}"}])
        nested.to_csv(OUT / "nested_selection.csv", index=False)
        pd.Series(chosen, name="model_chosen").to_csv(OUT / "nested_selection_choices.csv")
        log(nested.to_string(index=False))
        log("Models chosen per fold: " + pd.Series(chosen).value_counts().to_dict().__repr__())
    log(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
