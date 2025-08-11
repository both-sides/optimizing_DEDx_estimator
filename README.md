# ⚛️ Optimizing dE/dx Estimator — CMS, HSCP Searches 

**Purpose:**
This repository documents my work improving **dE/dx estimators** for HSCP (Heavy Stable Charged Particle) searches in CMS(Compact Muon Solenoid). It’s primarily a record of my methods, experiments, and results, rather than a plug-and-play package.

## 🔍 Project Overview

* Evaluating multiple per-track estimators of specific ionization (dE/dx) for a **1800 GeV gluino** sample (2018, 13 TeV CMS NanoAOD).
* Estimators benchmarked:

  * **Harmonic-2 mean**
  * **Adaptive truncated mean** (fraction removed scales with track length)
  * **Landau-MPV fit** using ROOT’s MINUIT/MIGRAD optimizer
* Focus areas:

  * Robustness to outliers and pathological hits
  * Improved reconstructed mass resolution and reduced bias
  * Diagnostic studies vs. track momentum (`p`) and pseudorapidity (`η`)

## 📈 Progress Highlights

* Landau-MPV fit converged successfully on **\~89.6%** of \~130k tracks.
* Observed a **secondary mass peak** in reconstructed mass distribution for Landau-MPV fit prototype — indicative of possible reliabilty after further optimization.
* Cumulative efficiency vs. dE/dx cut data saved to `cumulative_efficiencies.csv`.

## 🗂️ Repository Contents

```
optimizing_DEDx_estimator/
├─ notebooks/           # Estimator comparisons, mass reco, ROC curves
├─ src/                 # Estimator implementations, cleaning, fit utils
├─ output/              # Plots, tables, intermediate results
├─ env.yml              # Conda environment
├─ Makefile             # Automation tasks
└─ cumulative_efficiencies.csv
```

## 🛠 Tech Stack

* **ROOT / PyROOT**
* **Python**: NumPy, SciPy, Matplotlib
* **MINUIT / MIGRAD** (via ROOT)
* Jupyter notebooks for exploratory analysis
* Makefile for automation; Conda for environment management

## 📊 Current Methods

* **Harmonic-2 Mean**: robust to Landau tails compared to arithmetic mean.
* **Adaptive Truncated Mean**: % of high-ionizing hits dropped based on hit count.
* **Landau-MPV Fit**: per-track adaptively binned Landau fit, MPV used as dE/dx; convergence monitored and non-physical results filtered.
* **Data Cleaning**: isHighPurity selection, removal of Ih ≤ 0, filtering detector edge effects.

## 🗺 Roadmap

* [ ] Finalize systematic variation studies (truncation %, track selection)
* [ ] Publish clean estimator API in `src/`
* [ ] Complete ROC tables vs. CMS baselines
* [ ] Add bias and resolution summaries by kinematic bin

## 🔗 Context

* CMS HSCP searches
* PYTHIA v8.240, CP5 tune, NNPDF3.1
* √s = 13 TeV

*This README is a living document — progress and results will be updated as analysis continues.....*
