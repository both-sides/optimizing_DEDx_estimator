#!/usr/bin/env python3
import ROOT, csv
ROOT.gStyle.SetOptStat(0)
ROOT.EnableImplicitMT()   # optional speed-up on multicore
_KEEPALIVE = []


# -----------------------------
# Branch names (edit to match your ntuples)
# -----------------------------
HLT_BOOL         = "HLT_Mu50"
PT               = "IsoTrack_pt"
ETA              = "IsoTrack_eta"                        # if you don't have it, remove that cut below
FRAC_VALID       = "IsoTrack_fractionOfValidHits"
HIGH_PUR         = "IsoTrack_isHighPurityTrack"
NORM_CHI2        = "IsoTrack_normChi2"
DXY              = "IsoTrack_dxy"
DZ               = "IsoTrack_dz"
PIXELS_L2L4      = "DeDx_PixelNoL1NOM"
NDEDX            = None  # e.g., "IsoTrack_nDeDx"; set if available

# ---- Define the sequential *per-track* cuts (same track accumulates cuts) ----
# NOTE: Kept EXACTLY as in your working version (uses VecOps::abs)
track_cuts = [
    ("pT > 55 GeV",                   f"({PT} > 55)"),
    ("|eta| < 1",                     f"(ROOT::VecOps::abs({ETA}) < 1)"),
    ("# valid pixel hits L2-L4 >= 2", f"({PIXELS_L2L4} >= 2)"),
    ("Fraction of valid hits > 0.8",  f"({FRAC_VALID} > 0.8)"),
    # ("# dE/dx measurements >= 10",  f"({NDEDX} >= 10)") if NDEDX else None,
    ("High-purity track",             f"{HIGH_PUR}"),
    ("Track chi2/ndf < 5",            f"({NORM_CHI2} < 5)"),
    ("|dz| < 0.1 cm",                 f"(ROOT::VecOps::abs({DZ})  < 0.1)"),
    ("|dxy| < 0.02 cm",               f"(ROOT::VecOps::abs({DXY}) < 0.02)"),
]
track_cuts = [c for c in track_cuts if c is not None]

# -----------------------------
# Importable helper: build df_final WITHOUT argparse
# -----------------------------
def build_df_final(input_file: str, tree_path: str, weight_branch: str | None = None):
    """
    Minimal helper to obtain the final filtered RDataFrame without touching argparse.
    Returns (df_final, rows, cum_mask).
    """
    # --- open file / build RDataFrame (same logic as your script) ---
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {input_file}")
    obj = f.Get(tree_path)
    if obj and obj.InheritsFrom("TTree"):
        df = ROOT.RDataFrame(obj)      # TTree pointer
        _KEEPALIVE.append(f)           # <-- keep the TFile alive to avoid segfaults
    else:
        df = ROOT.RDataFrame(tree_path, input_file)

    # --- counting helper (weighted or unweighted) ---
    def count_events(dframe):
        return (float(dframe.Sum(weight_branch).GetValue())
                if weight_branch else int(dframe.Count().GetValue()))

    # All events
    N0 = count_events(df)
    rows = [("All events", N0, 1.0 if N0>0 else 0.0)]

    # Trigger only (event-level)
    dfe = df.Filter(HLT_BOOL, "Trigger")
    N1 = count_events(dfe)
    rows.append(("Trigger", N1, (N1/N0 if N0>0 else 0.0)))

    # Cumulative per-track mask
    cum_mask = None
    for label, expr in track_cuts:
        cum_mask = expr if cum_mask is None else f"({cum_mask}) && ({expr})"
        df_k = dfe.Filter(f"ROOT::VecOps::Any({cum_mask})", label)
        Nk = count_events(df_k)
        rows.append((label, Nk, (Nk/N0 if N0>0 else 0.0)))

    # Final selection DataFrame
    df_final = dfe.Filter(f"ROOT::VecOps::Any({cum_mask})", "Final selection")
    return df_final, rows, cum_mask


# -----------------------------
# CLI entrypoint (unchanged behavior when executed as a script)
# -----------------------------
def main():
    import argparse
    ap = argparse.ArgumentParser(description="Cumulative selection efficiency table")
    ap.add_argument("-i","--input", required=True, help="input ROOT file")
    ap.add_argument("-t","--tree",  required=True, help="tree name/path, e.g. HSCPFullAODAnalyzer/Events")
    ap.add_argument("--weight", default=None, help="event weight branch (optional)")
    ap.add_argument("--csv",    default="cumulative_efficiencies.csv", help="output CSV path")
    args = ap.parse_args()

    df_final, rows, cum_mask = build_df_final(args.input, args.tree, args.weight)

    # Print summary (same as before)
    def fmt_row(name, events, eff):
        return f"{name:<32} {events:>12}   {100.0*eff:6.3f} %"

    N_final = (df_final.Sum(args.weight).GetValue() if args.weight else df_final.Count().GetValue())
    print("\nEvents after full selection:", N_final)

    print("\nCUMULATIVE SELECTION EFFICIENCY\n")
    print(f"{'Selection':<32} {'Events':>12}   {'Eff. (cum to All)':>16}")
    print("-"*65)
    for name, events, eff in rows:
        print(fmt_row(name, events, eff))
    print("-"*65)

    with open(args.csv, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["Step", "CumulativeSelection", "Events", "Efficiency"])
        for i,(name, events, eff) in enumerate(rows):
            w.writerow([i, name, events, eff])
    print(f"Saved: {args.csv}")


if __name__ == "__main__":
    main()
