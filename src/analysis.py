import ROOT as rt
from ROOT import VecOps
import atexit
    
# open the ROOT file once, globally
root_file = rt.TFile.Open("../data/HSCPgluino_M-1800_fromAOD.root") 

# ensures it gets closed when the intergreter(jupyter kernel) closes
atexit.register(root_file.Close)


tree = root_file.Get("HSCPFullAODAnalyzer/Events")
if not tree:
    raise RuntimeError("Could not find HSCPFullAODAnalyzer/Events in the ROOT file")


# creates a dataframe globally
df = rt.RDataFrame(tree)

# filtered dataframe
df_filtered = df.Filter(
    "HLT_Mu50 && ROOT::VecOps::Sum(IsoTrack_pt > 55) > 0",
    "HLT_Mu50_and_anyIsoTrack55"
)


# some interested branch lists & color map ----------------------------------------------------------
HMNCSBR = ["DeDx_IhStrip1","DeDx_IhStrip", "DeDx_IhStrip3", "DeDx_IhStrip4"]
TRUNCSBR = ["DeDx_ItStrip0","DeDx_ItStrip5", "DeDx_ItStrip10", "DeDx_ItStrip15", "DeDx_ItStrip20", "DeDx_ItStrip25", "DeDx_ItStrip30", "DeDx_ItStrip35", "DeDx_ItStrip"]


COLOR_MAP = {
    1: rt.kViolet,   2: rt.kPink,   3: rt.kGray,   4: rt.kGreen,
    5: rt.kBlue,     6: rt.kYellow, 7: rt.kMagenta,8: rt.kCyan,
    9: rt.kOrange,  10: rt.kBlack
}


__all__ = ["root_file", "tree", "df", "df_filtered", "HscpData", "HMNCSBR", "TRUNCSBR", "COLOR_MAP"]


