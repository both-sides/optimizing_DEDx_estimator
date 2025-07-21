import math
import ROOT as rt
from typing import List, Dict, Tuple, Any, Optional
from datetime import date
import uuid
import numpy as np
from analysis import tree
import logging


## Constants ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TAG = date.today().isoformat()  # e.g. "2025-06-22"
OUTPUT_ROOT = "../output/root"
PLOTS_DIR   = "../output/plots"

COLOR_MAP = {
    1:  rt.kRed,
    2:  rt.kBlue,
    3:  rt.kGreen,
    4:  rt.kMagenta,
    5:  rt.kCyan,
    6:  rt.kOrange,
    7:  rt.kYellow,
    8:  rt.kViolet,
    9:  rt.kPink,
    10: rt.kAzure,
    11: rt.kSpring,
    12: rt.kTeal,
    13: rt.kGray,  
    14: rt.TColor.GetColor("#8B4513"),  # Brown
    15: rt.TColor.GetColor("#00CED1"),  # DarkTurquoise
}


FILTERS = [
    ("minPixelHits",   "DeDx_PixelNoL1NOM >= 2"),                    # Number of valid hits in L2-L4
    ("highPt",         "IsoTrack_pt > 55"),                          # pT > 55 GeV
    ("HLT_Mu50",       "HLT_Mu50"),                                  # Trigger requirement
    ("validHitsFrac",  "IsoTrack_fractionOfValidHits > 0.8"),       # Fraction of valid hits > 0.8
    ("minDedxHits",    "cluster_DeDxStrip.size() >= 10"),           # Number of dE/dx measurements >= 10
    ("highPurity",     "IsoTrack_isHighPurityTrack"),               # High purity track
    ("chi2",           "IsoTrack_normChi2 < 5"),                    # Normalized chi2 < 5
    ("dxy",            "abs(IsoTrack_dxy) < 0.02"),                 # |dxy| < 0.02 cm
    ("dz",             "abs(IsoTrack_dz) < 0.10"),                  # |dz| < 0.10 cm
]

## Functions ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
def apply_sequential_filters(rdf_in, filters=FILTERS, *, logger=None):
    """
    Apply filters one after another, printing (or logging) the surviving
    entry count after each cut. Returns the final RDataFrame.
    
    Parameters
    ----------
    rdf_in : ROOT.RDataFrame
        Input RDataFrame to filter
    filters : list of tuples
        List of (label, expression) pairs for filtering
    logger : logging.Logger, optional
        Logger to use for output
        
    Returns
    -------
    ROOT.RDataFrame
        Filtered RDataFrame
    """
    if logger is None:
        logger = logging.getLogger("FilterLog")
        if not logger.handlers:
            logging.basicConfig(level=logging.INFO,
                                format="%(message)s")

    counts = []            # Count() nodes – still lazy at this point
    rdf_curr = rdf_in

    # First count before any selections
    counts.append(("start", rdf_curr.Count()))
    
    logger.info(" ── Applying Sequential Filters ──")

    # Chain the filters with error handling
    for i, (label, expr) in enumerate(filters):
        try:
            logger.info(f"Applying filter {i+1}/{len(filters)}: {label}")
            logger.info(f"  Expression: {expr}")
            
            rdf_curr = rdf_curr.Filter(expr, label)
            counts.append((label, rdf_curr.Count()))
            
        except Exception as e:
            logger.error(f"  ❌ Filter '{label}' failed: {str(e)}")
            logger.error(f"     Expression: {expr}")
            
            # Try alternative expressions for common issues
            if "fabs" in expr:
                alt_expr = expr.replace("fabs", "abs")
                logger.info(f"  Trying alternative with abs(): {alt_expr}")
                try:
                    rdf_curr = rdf_curr.Filter(alt_expr, f"{label}_alt")
                    counts.append((f"{label}_alt", rdf_curr.Count()))
                    logger.info(f"  ✅ Alternative expression worked")
                except Exception as alt_e:
                    logger.error(f"  ❌ Alternative also failed: {str(alt_e)}")
                    logger.warning(f"  ⚠️  Skipping filter '{label}'")
                    continue
            else:
                logger.warning(f"  ⚠️  Skipping filter '{label}'")
                continue

    # Trigger the event loop exactly once
    logger.info(" ── Executing event loop ──")
    try:
        for label, counter in counts:
            n = counter.GetValue()          # evaluation happens here
            if label == "start":
                logger.info(f"{'Initial':<15s}: {n:10,d}")
            else:
                logger.info(f"{label:<15s}: {n:10,d}")

            if n == 0:
                logger.warning(f"   ↳ no entries left after '{label}' – "
                               "double-check this cut!")
                break
                
    except Exception as e:
        logger.error(f"Error during event loop execution: {str(e)}")
        logger.error("RDataFrame may be in an invalid state")
        raise

    return rdf_curr



def harmonic2_inloop(track):
    arr = np.asarray(track, dtype=np.float64)
    arr = arr[arr > 0]          # keep only positive hits
    if arr.size == 0:
        return 0.0              # or `np.nan`
    return np.sqrt(arr.size / np.sum(1.0 / arr**2))



def ArAvg_DEDx(cluster: List[List[float]]) -> List[List[float]]:
    """Arithmetic mean per track per event."""
    result = []
    for event in cluster:
        ev_avgs = []
        for track in event:
            if not track:
                ev_avgs.append(0.0)
            else:
                ev_avgs.append(sum(track) / len(track))
        result.append(ev_avgs)
    return result


def GeoAvg_DEDx(cluster: List[List[float]]) -> List[List[float]]:
    """Geometric mean per track per event."""
    result = []
    for event in cluster:
        ev_avgs = []
        for track in event:
            if not track:
                ev_avgs.append(0.0)
            else:
                ev_avgs.append(math.prod(track) ** (1.0 / len(track)))
        result.append(ev_avgs)
    return result


def h1Avg_DEDx(cluster: List[List[float]]) -> List[List[float]]:
    """Harmonic mean per track per event."""
    result = []
    for event in cluster:
        ev_avgs = []
        for track in event:
            if not track:
                ev_avgs.append(0.0)
            else:
                ev_avgs.append(len(track) / sum(1.0 / x for x in track))
        result.append(ev_avgs)
    return result


def records_equal(rec1, rec2, tol=1e-12):
    # compare every key/value in two dicts; tol for floating ≈ bit-identical
    if rec1.keys() != rec2.keys():
        return False
    for k in rec1:
        v1, v2 = rec1[k], rec2[k]
        if isinstance(v1, float):
            if abs(v1 - v2) > tol:
                return False
        else:
            if v1 != v2:
                return False
    return True



# builds a histogram stack and writes to the current open root file 
def write_stacked_histos(stack_name, hists, hists_title, canvas): # hists has to be a dictionary with {key = histogram name : value = histogram object (or pointers to that histogram object)}
    stack = rt.THStack(stack_name, hists_title)
    
    for hist in hists.values():
        stack.Add(hist) #pyroot stores histogram object pointers in the dictionary, I need to pull that out
    stack.Write()      #
    canvas.Write()
    return stack

# builds a histogram stack and writes to the current open root file 
def write_stacked_histos_ptr(stack_name, hists, hists_title, canvas): # hists has to be a dictionary with {key = histogram name : value = histogram object (or pointers to that histogram object)}
    stack = rt.THStack(stack_name, hists_title)
    
    for proxy in hists.values():
        stack.Add(proxy.GetPtr()) #pyroot stores histogram object pointers in the dictionary, I need to pull that out
    stack.Write()      #
    canvas.Write()
    return stack


#returns static start and end values adjusted with offset for each fit
def fit_range(lst: list, offset):
  start = min(lst) - offset
  end = max(lst) + offset
  return (start, end)


def seeds(hist):
  amp_guess  = hist.GetMaximum()                       # scale
  mpv_guess  = hist.GetBinCenter(hist.GetMaximumBin()) # MPV
  sigma_guess = 0.3 * hist.GetRMS() or 0.1*mpv_guess   # crude width
  
  return mpv_guess, amp_guess, sigma_guess


# Return (nbins, xmin, xmax) using the Freedman–Diaconis rule.
def freedman_diaconis_bins(values, *, range_pad=0.05):
    values = np.asarray(values, dtype=float)
    n = values.size
    if n < 2:
        raise ValueError("Need at least 2 points for IQR")
    q25, q75 = np.percentile(values, [25, 75])
    iqr = q75 - q25

    # Warn if IQR is zero
    if iqr == 0:
        print(f"[binning_warning] IQR=0 for {n} points "
              f"(min={values.min():.3g}, max={values.max():.3g})")


    h = 2.0 * iqr / n ** (1 / 3)          # Freedman–Diaconis width
    
    if h <= 0:                             # fallback for pathological IQR=0
        print(f"[binning_warning] Computed bin width h={h:.3g} ≤ 0; "
              "falling back to h=1e-12")
        h = 1e-12

    xmin, xmax = values.min(), values.max()
    span  = xmax - xmin

    nbins = int(np.ceil(span / h))
    if nbins < 1:
        print(f"[binning_error] nbins={nbins} < 1 (span={span:.3g}, h={h:.3g}); "
              "forcing nbins=1")
        nbins = 1

    # Padding edges a bit so min/max don't sit exactly on a bin edge
    pad = span * range_pad
    return max(nbins, 1), xmin - pad, xmax + pad


def fit(cluster: list, threshold: int = 2, *, max_hists: int = float('inf'), verbose=False):
    """
    Loop over tracks with len(track)>threshold, fit a Landau,
    and return dict of:
      - mpv_pt_pairs : list of (mpv, pt) tuples
      - params       : list of (mpv, sigma)
      - h2           : list of h2_mean
      - neg_mpvs     : list of "Event X Trk Y" with mpv<=0
    """
    
    # ensuring that the DEDx_IhStrip branch is on
    tree.SetBranchStatus("DeDx_IhStrip", 1)
    tree.SetBranchStatus("IsoTrack_pt", 1)
    
    # pre-allocate data containers
    mpv_pt_pairs = []     # (mpv, pt)
    params = []           # (mpv, sigma)
    harmonic2_means = []
    neg_ids = []
    neg_tracks = []
    
    # create one histogram & one TF1 and reuse them
    hist = rt.TH1F("h", "tmp", 1, 0, 1)     # binning/range will be reset
    f_landau = rt.TF1("f_landau", "landau", 0, 1)
    
    n_fits = 0
    logs = []  # buffer for fit printouts
    
    for event, tracks in zip(tree, cluster):
        for trk_idx, track in enumerate(tracks):
            # ---------- skip empties ----------
            if not track:
                continue
            
            if len(track) <= threshold:
                continue
            if n_fits >= max_hists:
                break
            
            # Get the pt value for this track
            try:
                pt_value = event.IsoTrack_pt[trk_idx]
            except (TypeError, IndexError):
                pt_value = event.IsoTrack_pt
            
            # ---------- single-hit ----------
            if len(track) == 1:
                mpv = track[0]
                sigma = 0.0
                harmonic2 = harmonic2_inloop(track)
            else:
                # using vectorised in loop calculation
                harmonic2 = harmonic2_inloop(track)
                
                # Freedman–Diaconis binning
                nbins, lo, hi = freedman_diaconis_bins(track)
                
                # reset & reconfigure the one histogram
                hist.Reset()
                hist.SetBins(nbins, 0, hi)
                bw = hist.GetBinWidth(nbins)
                hist.GetYaxis().SetTitle(f"Entries/{bw:.2f}")
                
                # filling the histogram
                for hit in track:
                    hist.Fill(hit)
                
                mpv_guess, amp_guess, sigma_guess = seeds(hist)
                
                f_landau.SetRange(0, hi)
                f_landau.SetParameters(amp_guess, mpv_guess, sigma_guess)  # seeding
                
                # I'm trying to keep Minuit away from crazy regions
                f_landau.SetParLimits(1, lo, hi)         # MPV must stay inside data
                f_landau.SetParLimits(2, 0.05, hi - lo)  # σ positive, < full range
                
                hist.Fit(f_landau, "RQ")

                # extracting results
                f = hist.GetFunction("f_landau")
                
                if not f:  # fit might have failed
                    continue
                
                mpv = f.GetParameter(1)
                sigma = f.GetParameter(2)
            
            n_fits += 1
            
            # Store MPV-pt pairs and other data
            mpv_pt_pairs.append((mpv, pt_value))
            params.append((mpv, sigma))
            harmonic2_means.append(harmonic2)
            
            if mpv <= 0:
                neg_ids.append(f"Event {event.event} Trk {trk_idx}")
                neg_tracks.append(track)
            
            if verbose:
                errs = [f.GetParError(i) for i in range(f.GetNpar())]
                logs.append(f" fit#{n_fits-1}: mpv = {mpv:.3g}±{errs[1]:.3g}  sigma = {sigma:.3g}±{errs[2]:.3g} pt = {pt_value:.3g}")
        
        if n_fits >= max_hists:
            break
    
    if verbose and logs:
        print("\n".join(logs))
    
    return {
        "mpv_pt_pairs": mpv_pt_pairs,
        "params": params,
        "h2": harmonic2_means,
        "neg_mpvs": neg_ids,
        "neg_tracks": neg_tracks
    }


def fit_mpv(cluster: list, threshold: int = 2, *, max_hists: int = float('inf'), verbose=False):
  """
  Loop over tracks with len(track)>threshold, fit a Landau,
  and return dict of:
    - corel   : list of (mpv, h2_mean)
    - params  : list of (mpv, sigma)
    - h2      : list of h2_mean
    - neg_mpvs: list of “Event X Trk Y” with mpv<=0
  """
  
  # ensuring that the DEDx_IhStrip branch is on
  tree.SetBranchStatus("DeDx_IhStrip", 1)
  
  # pre-allocate data containers
  corel_params   = []  # (mpv, h2_mean)
  params  = []  # (mpv, sigma)
  harmonic2_means= []
  neg_ids = []
  neg_tracks = []
  
  # create one histogram & one TF1 and reuse them
  hist     = rt.TH1F("h", "tmp",    1, 0, 1)     # binning/range will be reset
  f_landau= rt.TF1("f_landau", "landau", 0, 1)
  
  n_fits = 0
  logs   = []  # buffer for fit printouts
  
  for event, tracks in zip(tree, cluster):
    for trk_idx, track in enumerate(tracks): 
      # ---------- skip empties ----------
      if not track:
        continue
      
      if len(track) <= threshold:
        continue
      if n_fits >= max_hists:
        break
      
      # ---------- single-hit ----------
      if len(track) == 1:
        mpv = track[0]
        sigma = 0.0
        harmonic2 = harmonic2_inloop(track)
      else:
        #   #using precalculated track's h2 means
        #   try:
        #     harmonic2 = event.DeDx_IhStrip[trk_idx]
        #   except (TypeError, IndexError):
        #       # if it's a single float per track, omit the index
        #     harmonic2 = event.DeDx_IhStrip
        
        #using vectorised in loop calculation now
        harmonic2 = harmonic2_inloop(track)
        
        # Freedman–Diaconis binning
        nbins, lo, hi = freedman_diaconis_bins(track)
        
        # reset & reconfigure the one histogram
        hist.Reset()
        hist.SetBins(nbins, 0, hi)
        bw = hist.GetBinWidth(nbins)
        hist.GetYaxis().SetTitle(f"Entries/{bw:.2f}")
        
        # filling and drawing the hists
        for hit in track:
            hist.Fill(hit)
        
        mpv_guess, amp_guess, sigma_guess = seeds(hist)
        
        f_landau.SetRange(0, hi)
        f_landau.SetParameters(amp_guess, mpv_guess, sigma_guess) #seeding
        
        # I'm trying to keep Minuit away from crazy regions
        f_landau.SetParLimits(1, lo, hi)         # MPV must stay inside data
        f_landau.SetParLimits(2, 0.05, hi - lo)  # σ positive, < full range
        
        hist.Fit(f_landau, "RQ")

        # extracting results
        f = hist.GetFunction("f_landau")
        
        if not f:  # fit might have failed
          continue
        
        mpv = f.GetParameter(1)
        sigma = f.GetParameter(2)
      
      n_fits += 1
      
      corel_params.append((mpv, harmonic2))
      params.append((mpv, sigma))
      harmonic2_means.append(harmonic2)
      
      if mpv <= 0:
        neg_ids.append(f"Event {event.event} Trk {trk_idx}") 
        neg_tracks.append(track) # capture (the actual hits of these negative)tracks
  
      if verbose:
          errs = [f.GetParError(i) for i in range(f.GetNpar())]
          logs.append(f" fit#{n_fits-1}: mpv = {mpv:.3g}±{errs[1]:.3g}  sigma = {sigma:.3g}±{errs[2]:.3g} vs {harmonic2}")
   
        
    if n_fits >= max_hists:
      break
   
  if verbose and logs:
    print("\n".join(logs))
            
  return {
    "corel":    corel_params,
    "params":   params,
    "h2":       harmonic2_means,
    "neg_mpvs": neg_ids,
    "neg_tracks": neg_tracks
}

def draw_landau_fits(tree, cluster, threshold, max_hists):
    """
    Loop over events and tracks, build & fit Landau histograms, draw each fit,
    and return only the fit parameters, correlations, and the histograms themselves.

    Args:
      tree        : ROOT TTree or iterable of event objects
      clusters    : list-of-lists of hit‐value arrays, parallel to `tree`
      threshold   : minimum number of hits needed to attempt a fit
      max_hists   : maximum number of histograms/fits before stopping

    Returns:
      {
        "parameters": [(mpv, sigma), …],
        "corel":      [(mpv, harmonic2), …],
        "hists":      { hist_index: TH1F, … }
      }
    """
    parameters   = []
    corel_params = []
    hists        = {}
    hist_count   = 0

    # single canvas reused for speed
    canvas = rt.TCanvas("c_landau", "Landau Fits", 800, 600)

    for event, tracks in zip(tree, cluster):
        for trk_idx, track in enumerate(tracks):
            if len(track) <= threshold or hist_count >= max_hists:
                continue

            # get the precomputed harmonic‐2 mean from the event
            try:
                h2 = event.DeDx_IhStrip[trk_idx]
            except (TypeError, IndexError):
                h2 = event.DeDx_IhStrip
            
            # #using vectorised in loop calculation now
            # h2 = harmonic2_inloop(track)

            # binning by Freedman‐Diaconis
            nbins, lo, hi = freedman_diaconis_bins(track)

            # fill histogram
            name  = f"h{hist_count}"
            title = f"Event {event.event}, Trk {trk_idx};dE/dx (MeV/cm);Entries"
            hist  = rt.TH1F(name, title, nbins, 0, hi)
            width = hist.GetBinWidth(1)
            hist.GetYaxis().SetTitle(f"Entries/{width:.2f}")

            for hit in track:
                hist.Fill(hit)

            # seed the Landau
            mpv0, amp0, sigma0 = seeds(hist)
            f_landau = rt.TF1("f_landau", "landau", 0, hi)
            f_landau.SetParameters(amp0, mpv0, sigma0)
            f_landau.SetParLimits(1, lo, hi)         # MPV inside data
            f_landau.SetParLimits(2, 0.01, hi - lo)  # σ > 0

            # do the fit quietly and store the hist
            hist.Fit(f_landau, "RQ")
            hists[hist_count] = hist

            # draw result
            canvas.cd()
            hist.Draw("P")
            f_landau.Draw("same")
            canvas.Update()

            # pull out the fit results
            mpv   = f_landau.GetParameter(1)
            sigma = f_landau.GetParameter(2)
            parameters.append((mpv, sigma))
            corel_params.append((mpv, h2))

            hist_count += 1
            if hist_count >= max_hists:
                break
        if hist_count >= max_hists:
            break

    return {
        "parameters": parameters,
        "corel":      corel_params,
        "hists":      hists,
    }



def build_event_index(tree: rt.TTree, event_branch: str = "event") -> Dict[int, int]:
    """
    Scan `tree` once and return a dict mapping
      event_number (int) -> tree entry index (0..N-1)
    """
    idx = {}
    for eidx in range(tree.GetEntries()):
        tree.GetEntry(eidx)
        evt = int(getattr(tree, event_branch))  
        idx[evt] = eidx
    return idx


def get_attrs_for_labels(tree: rt.TTree,
                         labels: List[str], 
                         attrs: List[str], 
                         event_to_entry: Optional[Dict[int, int]] = None,
                         event_branch: str = "event")-> Tuple[Dict[str, Dict[str, Any]], Dict[int, int]]:
    
    """
    Parameters
    ----------
    tree            : ROOT TTree
    labels          : ["Event <num> Trk <idx>", ...]
    attrs           : list of branch names you want (e.g. ['DeDx_IhStrip','Isotrack_pt'])
    event_to_entry  : cache from build_event_index(); if None we build it
    event_branch    : name of the event-number branch (default 'event')

    Returns
    -------
    results         : { label : { attr : value, ... }, ... }
    event_to_entry  : the (possibly newly built) cache you can reuse
    """
 
    # build cache on first call
    if event_to_entry is None:
        event_to_entry = build_event_index(tree, event_branch)

    results: Dict[str, Dict[str, Any]] = {}
    for label in labels:
        # parse out the numbers
         # "Event 43249 Trk 2"
        try:
            _, evt_str, _, trk_str = label.split()
            evt_num  = int(evt_str)
            trk_idx  = int(trk_str)
        except ValueError:
            raise ValueError(f"Label '{label}' not in 'Event <num> Trk <idx>' format")

        entry = event_to_entry.get(evt_num)
        if entry is None:
            results[label] = {a: None for a in attrs}
            continue

        # load that one event
        tree.GetEntry(entry)

        per_track_data = {}
        for a in attrs:
            raw = getattr(tree, a)
            # try vector-style access, else scalar
            try:
                per_track_data[a] = raw[trk_idx]
            except (TypeError, IndexError):
                per_track_data[a] = raw
        results[label] = per_track_data

    return results, event_to_entry
    

def report_integrals(hists, mpv_hist=None):
    """
    Print total entries (incl. under/overflow) for each track histogram,
    and for the MPV histogram if provided.
    
    Args:
      hists    : dict of { idx : TH1F }
      mpv_hist : TH1F of MPV values (optional)
    """
    for idx, hist in hists.items():
        nb = hist.GetNbinsX()
        total = hist.Integral(0, nb+1)
        print(f"[hist #{idx}] bins=0→{nb+1}  entries={total:g}")

    if mpv_hist is not None:
        nb = mpv_hist.GetNbinsX()
        total = mpv_hist.Integral(0, nb+1)
        print(f"[MPV hist] bins=0→{nb+1}  entries={total:g}")


def check_branch_shapes(data_dict, keys, *, ref_key=None, verbose=True):
    """
    Verify that the given branches have matching shapes row by row.

    Parameters
    ----------
    data_dict : dict
        Output of RDataFrame.AsNumpy()  maps column name → NumPy array.
    keys : Sequence[str]
        Branch names to compare.
    ref_key : str, optional
        Reference branch for per row inner length comparison.
        Defaults to the first key in `keys`.
    verbose : bool, default True
        If True, print a summary and first mismatched rows lengths.

    Returns
    -------
    ok : bool
        True if every row is consistent across all branches.
    bad_idx : list[int]
        Indices of rows with mismatched inner lengths (empty when ok=True).
    """
    if not keys:
        raise ValueError("`keys` cannot be empty.")
    if ref_key is None:
        ref_key = keys[0]

    # ---------- top‑level length check ----------
    n_rows = len(data_dict[ref_key])
    for k in keys:
        if len(data_dict[k]) != n_rows:
            raise ValueError(
                f"Top level length mismatch for {k}: {len(data_dict[k])} vs {n_rows}"
            )

    # ---------- row‑by‑row inner‑length check ----------
    row_ok = [
        all(len(data_dict[k][i]) == len(data_dict[ref_key][i]) for k in keys)
        for i in range(n_rows)
    ]
    bad_idx = [i for i, ok in enumerate(row_ok) if not ok]

    if verbose:
        print(f"All rows consistent? {len(bad_idx) == 0}")
        if bad_idx:
            print(f"Found {len(bad_idx)} mismatched rows: {bad_idx[:10]} …")
            # Inspect the first problematic row
            i = bad_idx[0]
            for k in keys:
                print(f"{k}: len = {len(data_dict[k][i])}")

    return len(bad_idx) == 0, bad_idx



def flatten_branches(data_dict, branches, *,
                     drop_sentinel=None,
                     return_stats=True,
                     verbose=True):
    """
    Flatten variable–length RDataFrame branches that were pulled in with
    `AsNumpy`, optionally drop a sentinel value (e.g. -1), and (optionally)
    return or print simple stats.

    Parameters
    ----------
    data_dict : dict
        Output of RDataFrame.AsNumpy() – column name → NumPy array (dtype=object
        for std::vector branches, numeric dtype for scalar branches).

    branches : str | Sequence[str]
        One branch name or a list/tuple of names to flatten.

    drop_sentinel : float | int | None, default None
        Value that marks invalid entries (e.g. -1.0 in DeDx vectors).
        If given, those elements are filtered out after flattening.

    return_stats : bool, default True
        If True, returns a dict keyed by branch with (min, max, mean, n)
        so you can log or assert ranges in calling code.

    verbose : bool, default True
        Print shapes and basic stats to stdout.

    Returns
    -------
    flat_dict : dict[str, np.ndarray]
        Branch → 1‑D NumPy array of flattened (and cleaned) values.

    stats_dict : dict[str, tuple]  (only if return_stats=True)
        Branch → (min, max, mean, length)
    """
    if isinstance(branches, str):
        branches = [branches]

    flat_dict  = {}
    stats_dict = {}

    for br in branches:
        col = data_dict[br]

        # Flatten: concatenate list/array elements for each event
        if col.dtype == "O":        # vector branch
            flat = np.concatenate(
                [np.asarray(ev) for ev in col if len(ev) > 0]
            )
        else:                       # already scalar array
            flat = col

        # Drop sentinel values, if requested
        if drop_sentinel is not None:
            flat = flat[flat != drop_sentinel]

        flat_dict[br] = flat

        if return_stats:
            stats = (float(flat.min()), float(flat.max()),
                     float(flat.mean()), len(flat))
            stats_dict[br] = stats

        if verbose:
            msg = (f"{br}: shape={flat.shape}, "
                   f"range={flat.min():.3f}–{flat.max():.3f}")
            if drop_sentinel is not None:
                msg += f"   (sentinel {drop_sentinel} dropped)"
            print(msg)

    return (flat_dict, stats_dict) if return_stats else flat_dict



## Classes ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#TH1s histogram drawer class
class HistogramDrawer:
    def __init__(self,
                 canvas_width: int = 800,
                 canvas_height: int = 600):
                #  dpi: int = 100):
        #initialize default drawing params
        self.canvas_w = canvas_width
        self.canvas_h = canvas_height
        #self.dpi = dpi
        
    #create and return a fresh TCanvas   
    def _new_canvas(self, 
                    name: str,
                    title: str,
                    width: int = None,
                    height: int = None) -> rt.TCanvas:
       
        w = width or self.canvas_w
        h = height or self.canvas_h
        
        #close existing canvas of same name so JSROOT will re-render
        existing = rt.gROOT.FindObject(name)
        if existing:
            existing.Close()
        
        return rt.TCanvas(name, title, w, h)
    #set axis titles    
    def set_axis_titles(self,
                        hist: rt.TH1,
                        xlabel: str = "",
                        ylabel: str = "Number of tracks"):
        hist.GetXaxis().SetTitle(xlabel)
        hist.GetYaxis().SetTitle(ylabel)
    
    #Draw a single histogram on its own canvas    
    def draw_hist(self,
                  hist: rt.TH1,
                  canvas_name: str = "canvas",
                  canvas_title: str = "",
                  xlabel: str = "",
                  ylabel: str = "Number of tracks",
                  logy: bool = False,
                  draw_opts: str = "hist"):
        
        c = self._new_canvas(canvas_name, canvas_title)
        if logy:
            c.SetLogy(1)
        self.set_axis_titles(hist, xlabel, ylabel)
        hist.Draw(draw_opts)
        c.Update()
        return c
    
    #create a THStack and draw it
    def draw_stack(self,
                   hists: list[rt.TH1],
                   stack_name: str = "stack",
                   canvas_title: str = "",
                   xlabel: str = "",
                   ylabel: str = "Number of tracks",
                   logy: bool = False,
                   draw_opts: str = ""):
        
        stack = rt.THStack(stack_name, canvas_title)
        for i, hist in enumerate(hists):
            hist.SetLineColor(COLOR_MAP[i+1])
            stack.Add(hist)
        
        c = self._new_canvas(f"{stack_name}'s canvas", canvas_title)
        if logy:
            c.SetLogy(1)
            
        stack.GetXaxis().SetTitle(xlabel)
        stack.GetYaxis().SetTitle(ylabel)
        
        stack.Draw(draw_opts)
        c.Update()
        return (c, stack)
    
    
    #TODO: THE LOWER PAD IS SOMEHOW NOT RENDERING CORECTLY
    
    #Two‐panel ratio plot: top = overlay, bottom = ratio.
    # def draw_ratio(self,
    #                h_num: rt.TH1,
    #                h_den: rt.TH1,
    #                canvas_name: str = "ratio",
    #                title: str = "",
    #                xlabel: str = "",
    #                ylabel: str = "Number of tracks",
    #                ratio_ylabel: str = "Ratio",
    #                ratio_range: tuple[float, float] = (-5, 5),
    #                logy: bool = False) -> rt.TCanvas:
    #     """
    #     Draws a comparison of two histograms (h_num over h_den) with a ratio pad below.
    #     Returns the ROOT TCanvas containing the two pads.
    #     """
    #     # Close any existing canvases so ROOT will make a new one
    #     # for canv in rt.gROOT.GetListOfCanvases():
    #     #     canv.Close()

    #     # 1) Create canvas
    #     c = self._new_canvas(canvas_name, title or "Ratio Plot", 800, 800)

    #     # 2) Create pads: top (0.3–1) and bottom (0–0.3)
    #     pad1 = rt.TPad(f"{canvas_name}_pad1", "Top pad", 0, 0.3, 1, 1)
    #     pad2 = rt.TPad(f"{canvas_name}_pad2", "Bottom pad", 0, 0, 1, 0.3)
    #     pad1.SetBottomMargin(0)
    #     pad2.SetTopMargin(0.02)
    #     pad2.SetBottomMargin(0)
    #     pad1.Draw()
    #     pad2.Draw()

    #     # 3) Draw histograms on top pad
    #     pad1.cd()
    #     if logy:
    #         pad1.SetLogy(1)

    #     h_den.SetLineColor(rt.kBlack)
    #     h_num.SetLineColor(rt.kRed)
    #     maxval = max(h_den.GetMaximum(), h_num.GetMaximum())
    #     h_den.SetMaximum(maxval * 1.2)

    #     h_den.SetTitle(title)
    #     h_den.GetYaxis().SetTitle(ylabel)
    #     h_den.Draw()
    #     h_num.Draw("same")
    #     pad1.Update()

    #     # 4) Build ratio in bottom pad
    #     pad2.cd()
    #     ratio = h_num.Clone(f"{h_num.GetName()}_ratio")
    #     ratio.Divide(h_den)
    #     ratio.SetMarkerStyle(20)
    #     ratio.SetTitle("")
    #     ratio.GetYaxis().SetTitle(ratio_ylabel)
    #     ratio.GetXaxis().SetTitle(xlabel)

    #     # Adjust label sizes
    #     ratio.GetYaxis().SetTitleSize(0.05)
    #     ratio.GetYaxis().SetTitleOffset(0.4)
    #     ratio.GetYaxis().SetLabelSize(0.05)
    #     ratio.GetXaxis().SetTitleSize(0.08)
    #     ratio.GetXaxis().SetLabelSize(0.08)

    #     ratio.SetMinimum(ratio_range[0])
    #     ratio.SetMaximum(ratio_range[1])
    #     ratio.Draw("EP")
    #     pad2.Update()

    #     # 5) Finalize canvas
    #     c.cd()
    #     c.Modified()
    #     c.Draw()

    #     return c
    
    def create_two_ratio(self,
                       h_num: rt.TH1,
                       h_den: rt.TH1,
                       option: str = "pois",
                       title: str = "") -> rt.TRatioPlot:
        """
        Build and return a TRatioPlot for h_num / h_den (no drawing).
        """
        
        #creating a clone of h_num..... i'd be formatting some of its attributes
        h_num = h_num.Clone(f"{h_num.GetName()} vs {h_den.GetName()}")
        h_num.SetTitle(title)
        h_den.SetLineColor(rt.kBlack)
        h_num.SetLineColor(rt.kRed)
        rp = rt.TRatioPlot(h_num, h_den, option)
        return rp

    def draw_ratio(self,
                   rp: rt.TRatioPlot,
                   canvas_name: str = "c",
                   title: str = "",
                   xlabel: str = "DEDx [MeV]",
                   ylabel: str = "Number of tracks",
                   yratio_label: str = "h_num/h_den",
                   range: list = [-5, 5],
                   logy: bool = False) -> rt.TCanvas:
 
        # Close any existing canvas with this name
        for canv in list(rt.gROOT.GetListOfCanvases()):
            if canv.GetName() == canvas_name:
                canv.Close()

        # Create new canvas
        c = self._new_canvas(canvas_name, title)
        if logy:
            c.SetLogy()
            
        rp.SetH1DrawOpt("hist")      # draw first histogram as an empty‐filled outline
        rp.SetH2DrawOpt("hist same") # draw second histogram with same binning
        rp.SetGraphDrawOpt("P")

        # Draw ratio plot
        rp.Draw()
        
        c.Update()
        # — Upper pad X–axis:
        ux = rp.GetUpperRefXaxis()
        ux.SetTitleOffset(1.0)
        ux.SetLabelSize(0.04)
        
        
        # lower x axis
        lx = rp.GetLowerRefXaxis()
        lx.SetTitle(xlabel)
        
        # — Upper pad Y–axis:
        uy = rp.GetUpperRefYaxis()
        uy.SetTitle(ylabel)    # the histogram count label
        uy.SetTitleOffset(1.2)
        uy.SetLabelSize(0.04)
        
        # — Lower pad Y–axis (the ratio):
        ly = rp.GetLowerRefYaxis()
        ly.SetTitle(yratio_label)             
        ly.SetTitleOffset(0.8)
        ly.SetLabelSize(0.02)

        lg = rp.GetLowerRefGraph()
        lg.SetMinimum(range[0])
        lg.SetMaximum(range[1])
        rp.GetLowerPad().Modified()
        rp.GetLowerPad().Update() 
        
        c.Modified()
        c.Update()
        
        return c


        
        


#TODO: include the  fit residual method(ratio of a fit and the histogram), and THStack vs TH1 ratio plot method


    def save(self,
             canvas: rt.TCanvas,
             filename: str,
             formats: list[str] = ("png", "pdf")):
        """
        Save a canvas to disk in one or more formats.
        """
        for fmt in formats:
            canvas.SaveAs(f"{filename}.{fmt}")








