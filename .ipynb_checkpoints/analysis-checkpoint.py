import ROOT as rt
from ROOT import VecOps
import numpy as np
import matplotlib.pyplot as plt
import math
import array

f = rt.TFile.Open("data/HSCPgluino_M-1800_fromAOD.root")
dir = f.GetDirectory("HSCPFullAODAnalyzer")

obj = dir.Get("Events")
#print(obj.ClassName())   

tree = f.Get("HSCPFullAODAnalyzer/Events")

#evnt_cut_flow_func = dir.Get("EventCutFlow")
#hscp_cut_flow_func = dir.Get("HSCPCutFlow")

filename = "data/HSCPgluino_M-1800_fromAOD.root"
treename = "HSCPFullAODAnalyzer/" + dir.Get("Events").GetName()


df0 = rt.RDataFrame(treename, filename)

df = df0.Filter(
    "HLT_Mu50 && ROOT::VecOps::Sum(IsoTrack_pt > 55) > 0",
    "HLT_Mu50_and_anyIsoTrack55"
)

# dictionaries to hold some interested branch names and created histogram names
hmncSbr = ["DeDx_IhStrip1","DeDx_IhStrip", "DeDx_IhStrip3", "DeDx_IhStrip4"]
hists1 = {}
truncSbr = ["DeDx_ItStrip0","DeDx_ItStrip5", "DeDx_ItStrip10", "DeDx_ItStrip15", "DeDx_ItStrip20", "DeDx_ItStrip25", "DeDx_ItStrip30", "DeDx_ItStrip35", "DeDx_ItStrip"]
hists2 = {}


color_map = {
    1: rt.kViolet,
    2: rt.kPink,
    3: rt.kGray,
    4: rt.kGreen,
    5: rt.kBlue,
    6: rt.kYellow,
    7: rt.kMagenta,
    8: rt.kCyan,
    9: rt.kOrange,
    10: rt.kBlack
}

