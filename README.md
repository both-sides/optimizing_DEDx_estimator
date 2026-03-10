# Optimizing dE/dx Estimators — CMS HSCP Searches

This repository documents work on improving dE/dx estimators for Heavy Stable Charged Particle (HSCP) searches in CMS. It is primarily a record of methods, experiments, and results — not a plug-and-play package.

---

## Overview

Evaluating per-track estimators of specific ionization (dE/dx) on a 1800 GeV gluino sample (2018, 13 TeV CMS NanoAOD).

**Estimators benchmarked:**
- Harmonic-2 mean
- Adaptive truncated mean (fraction removed scales with track length)
- Landau-MPV fit via ROOT's MINUIT/MIGRAD optimizer

**Focus areas:**
- Robustness to outliers and pathological hits
- Improved reconstructed mass resolution and reduced bias
- Diagnostic studies vs. track momentum (`p`) and pseudorapidity (`η`)

---

## Results

- Landau-MPV fit converged on ~89.6% of ~130k tracks
- Observed a secondary mass peak in the reconstructed mass distribution for the Landau-MPV prototype — indicative of possible reliability gains after further optimization
- Cumulative efficiency vs. dE/dx cut data saved to `cumulative_efficiencies.csv`

---

## Repository Structure
```
optimizing_DEDx_estimator/
├─ notebooks/           # Estimator comparisons, mass reco, ROC curves
├─ src/                 # Estimator implementations, cleaning, fit utils
├─ output/              # Plots, tables, intermediate results
├─ data/
├─ env.yml              # Conda environment
├─ Makefile             # Automation tasks
└─ cumulative_efficiencies.csv
```

---

## Methods

**Harmonic-2 Mean** — robust to Landau tails compared to arithmetic mean.

**Adaptive Truncated Mean** — fraction of high-ionizing hits dropped scales with hit count per track.

**Landau-MPV Fit** — per-track adaptively binned Landau fit; MPV used as dE/dx estimate. Convergence monitored; non-physical results filtered.

**Data Cleaning** — `isHighPurity` selection, removal of Ih ≤ 0, filtering of detector edge effects.

---

## Stack

- ROOT / PyROOT, MINUIT/MIGRAD
- Python: NumPy, SciPy, Matplotlib
- Jupyter notebooks, Makefile, Conda

---

## Roadmap

- [ ] Finalize systematic variation studies (truncation %, track selection)
- [ ] Publish clean estimator API in `src/`
- [ ] Complete ROC tables vs. CMS baselines
- [ ] Add bias and resolution summaries by kinematic bin

---

## Context

CMS HSCP searches · PYTHIA v8.240, CP5 tune, NNPDF3.1 · √s = 13 TeV

*Living document — updated as analysis continues.*
