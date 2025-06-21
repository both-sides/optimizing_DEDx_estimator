import ROOT as rt
from ROOT import VecOps
import numpy as np
import matplotlib.pyplot as plt
import math
import array

root_file = rt.TFile.Open("../data/HSCPgluino_M-1800_fromAOD.root")
hs_analyzer_dir = root_file.GetDirectory("HSCPFullAODAnalyzer")

events_obj = hs_analyzer_dir.Get("Events")
#print(obj.ClassName())   

tree = root_file.Get("HSCPFullAODAnalyzer/Events")

#evnt_cut_flow_func = hs_analyzer_dir.Get("EventCutFlow")
#hscp_cut_flow_func = hs_analyzer_dir.Get("HSCPCutFlow")

filename = "../data/HSCPgluino_M-1800_fromAOD.root"
treename = "HSCPFullAODAnalyzer/" + events_obj.GetName()


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

#color map
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

#Arithmetic everage for each DEDx for each track
def ArAvg_DEDx(cluster):
    Aavg_DEDx = []
    ## each event should be an array of averages where each element in that array corresponds to the
    ## respective track average at that index
    zero_count = 0
    for i, event_tracks in enumerate(cluster):
        # tracks also coresponds to each event
        #each track's average is a value
        tracks_ar_avg = [] # array to hold all track DEDx averages in each tracks(event)
        for track in event_tracks:
            if len(track) == 0:
                track_avg = 0
                tracks_ar_avg.append(track_avg)
                zero_count += 1
            else:
                track_avg = sum(track)/len(track)
                tracks_ar_avg.append(track_avg)
        Aavg_DEDx.append(tracks_ar_avg)
    return Aavg_DEDx


#Geometric everage for each DEDx for each track
def GeoAvg_DEDx(cluster):
    Geo_avg_DEDx = []
    ## each event should be an array of averages where each element in that array corresponds to the
    ## respective track average at that index
    zero_count = 0
    for i, event_tracks in enumerate(cluster):
        # tracks also coresponds to each event
        #each track's average is a value
        tracks_geo_avg = [] # array to hold all track DEDx averages in each tracks(event)
        for track in event_tracks:
            if len(track) == 0: # no hits >> define geo-mean = 0
                tracks_geo_avg.append(0)
                zero_count += 1
            else:
                #if any DEDx == 0, product = 0, geo-mean = 0 >> should that make sense here
                track_avg = math.prod(track)**(1/len(track)) 
                tracks_geo_avg.append(track_avg)
        Geo_avg_DEDx.append(tracks_geo_avg)
    return Geo_avg_DEDx


#Harmonic-1 everage for each DEDx for each track
def h1Avg_DEDx(cluster):
    h1_avg_DEDx = []
    ## each event should be an array of averages where each element in that array corresponds to the
    ## respective track average at that index
    zero_count = 0
    for i, event_tracks in enumerate(cluster):
        # tracks also coresponds to each event
        #each track's average is a value
        tracks_h1_avg = [] # array to hold all track DEDx averages in each tracks(event)
        for track in event_tracks:
            if len(track) == 0: #or any(x == 0 for x in track) (seems that there's no dEdx value that is 0): # no hits >> define harmonic mean = 0
                tracks_h1_avg.append(0)
                zero_count += 1
            else:   
                track_avg = (len(track))/sum(x**(-1) for x in track)
                tracks_h1_avg.append(track_avg)
        h1_avg_DEDx.append(tracks_h1_avg)
    return h1_avg_DEDx
