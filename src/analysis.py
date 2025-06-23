import ROOT as rt
from ROOT import VecOps
import atexit
# import numpy as np
# import matplotlib.pyplot as plt
# import math
# import array

#-------------------------------------------------------------------------------------------
#  context manager for input file & TTree
# class HscpData:
#     """Open the AOD file and expose .tree Closes on exit"""
#     def __init__(self,
#                  path: str = "../data/HSCPgluino_M-1800_fromAOD.root",
#                  treename: str = "HSCPFullAODAnalyzer/Events"):
#         # self.file = rt.TFile.Open(path)
#         # self.tree = self.file.Get("HSCPFullAODAnalyzer/Events")
#         self._path     = path
#         self._treename = treename
#         self.file = None
#         self.tree = None
        
#     def __enter__(self):
#         self.file = rt.TFile.Open(self._path)
#         if not self.file or not self.file.IsOpen():
#             raise IOError(f"Cannot open {self._path}")
#         self.tree = self.file.Get(self._treename)
#         if not self.tree:
#             raise KeyError(f"Tree {self._treename} not found")
#         return self
    
#     # def close(self):
#     #     if self.file and self.file.IsOpen():
#     #         self.file.Close()
    
#     def __exit__(self, exc_type, exc_val, exc_tb):
#         if self.file and self.file.IsOpen():
#             self.file.Close()
#     # # context-manager hooks
#     # def __enter__(self): return self
#     # def __exit__(self, exc_type, exc, tb): self.close()
    
# open the ROOT file once, globally
root_file = rt.TFile.Open("../data/HSCPgluino_M-1800_fromAOD.root") 

# ensures it gets closed when the intergreter(jupyter kernel) closes
atexit.register(root_file.Close)

# hs_analyzer_dir = root_file.GetDirectory("HSCPFullAODAnalyzer")

# events_obj = hs_analyzer_dir.Get("Events")
# #print(obj.ClassName())   

tree = root_file.Get("HSCPFullAODAnalyzer/Events")
if not tree:
    raise RuntimeError("Could not find HSCPFullAODAnalyzer/Events in the ROOT file")

# #evnt_cut_flow_func = hs_analyzer_dir.Get("EventCutFlow")
# #hscp_cut_flow_func = hs_analyzer_dir.Get("HSCPCutFlow")

# filename = "../data/HSCPgluino_M-1800_fromAOD.root"
# treename = "HSCPFullAODAnalyzer/" + events_obj.GetName()

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


