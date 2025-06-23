import math
import ROOT as rt
from typing import List
from datetime import date


## Constants ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TAG = date.today().isoformat()  # e.g. "2025-06-22"
OUTPUT_ROOT = "../output/root"
PLOTS_DIR   = "../output/plots"


## Functions ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# builds a histogram stack and writes to the current open root file 
def write_stacked_histos(stack_name, hists, hists_title, canvas): # hists has to be a dictionary with {key = histogram name : value = histogram object (or pointers to that histogram object)}
    stack = rt.THStack(stack_name, hists_title)
    
    for proxy in hists.values():
        stack.Add(proxy.GetPtr()) #pyroot stores histogram object pointers in the dictionary, I need to pull that out
    stack.Write()      #
    canvas.Write()
    