import ROOT as rt
from ROOT import VecOps
import atexit
    
# open the ROOT file once, globally
root_file = rt.TFile.Open("../data/HSCPgluino_M-1800_fromAOD.root") 
data_root_file = rt.TFile.Open("../data/Data18_03639171301F_fromAOD.root")
# ensures it gets closed when the intergreter(jupyter kernel) closes
atexit.register(root_file.Close)


tree = root_file.Get("HSCPFullAODAnalyzer/Events")
if not tree:
    raise RuntimeError("Could not find HSCPFullAODAnalyzer/Events in the ROOT file")

tree1 = data_root_file.Get("HSCPFullAODAnalyzer/Events")
if not tree1:
    raise RuntimeError("Could not find HSCPFullAODAnalyzer/Events in the ROOT file")


# creates a dataframe globally
df = rt.RDataFrame(tree)
df1 = rt.RDataFrame(tree1)

sel = """
  HLT_Mu50
  && ROOT::VecOps::Any(IsoTrack_pt                     > 55)
  && ROOT::VecOps::Any(IsoTrack_fractionOfValidHits   > 0.8)
  && ROOT::VecOps::Any(IsoTrack_isHighPurityTrack)
  && ROOT::VecOps::Any(IsoTrack_normChi2              < 5)
  && ROOT::VecOps::Any(abs(IsoTrack_dxy)              < 0.02)
  && ROOT::VecOps::Any(abs(IsoTrack_dz)               < 0.1)
  && ROOT::VecOps::Any(DeDx_PixelNoL1NOM              >= 2)
  
"""

# filtered dataframe
df_filtered = df.Filter(
    sel, "all cuts"
)


# some interested branch lists & color map ----------------------------------------------------------
HMNCSBR = ["DeDx_IhStrip1","DeDx_IhStrip", "DeDx_IhStrip3", "DeDx_IhStrip4"]
TRUNCSBR = ["DeDx_ItStrip0","DeDx_ItStrip5", "DeDx_ItStrip10", "DeDx_ItStrip15", "DeDx_ItStrip20", "DeDx_ItStrip25", "DeDx_ItStrip30", "DeDx_ItStrip35", "DeDx_ItStrip"]


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
    13: rt.kGray,  # neutral mid‐tone
    # two extra custom colors:
    14: rt.TColor.GetColor("#8B4513"),  # Brown
    15: rt.TColor.GetColor("#00CED1"),  # DarkTurquoise
}


__all__ = ["root_file", "tree", "df", "df_filtered", "HscpData", "HMNCSBR", "TRUNCSBR", "COLOR_MAP"]


