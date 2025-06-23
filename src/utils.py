import math
from typing import List


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








# #Arithmetic everage for each DEDx for each track
# def ArAvg_DEDx(cluster):
#     Aavg_DEDx = []
#     ## each event should be an array of averages where each element in that array corresponds to the
#     ## respective track average at that index
#     zero_count = 0
#     for i, event_tracks in enumerate(cluster):
#         # tracks also coresponds to each event
#         #each track's average is a value
#         tracks_ar_avg = [] # array to hold all track DEDx averages in each tracks(event)
#         for track in event_tracks:
#             if len(track) == 0:
#                 track_avg = 0
#                 tracks_ar_avg.append(track_avg)
#                 zero_count += 1
#             else:
#                 track_avg = sum(track)/len(track)
#                 tracks_ar_avg.append(track_avg)
#         Aavg_DEDx.append(tracks_ar_avg)
#     return Aavg_DEDx


# #Geometric everage for each DEDx for each track
# def GeoAvg_DEDx(cluster):
#     Geo_avg_DEDx = []
#     ## each event should be an array of averages where each element in that array corresponds to the
#     ## respective track average at that index
#     zero_count = 0
#     for i, event_tracks in enumerate(cluster):
#         # tracks also coresponds to each event
#         #each track's average is a value
#         tracks_geo_avg = [] # array to hold all track DEDx averages in each tracks(event)
#         for track in event_tracks:
#             if len(track) == 0: # no hits >> define geo-mean = 0
#                 tracks_geo_avg.append(0)
#                 zero_count += 1
#             else:
#                 #if any DEDx == 0, product = 0, geo-mean = 0 >> should that make sense here
#                 track_avg = math.prod(track)**(1/len(track)) 
#                 tracks_geo_avg.append(track_avg)
#         Geo_avg_DEDx.append(tracks_geo_avg)
#     return Geo_avg_DEDx


# #Harmonic-1 everage for each DEDx for each track
# def h1Avg_DEDx(cluster):
#     h1_avg_DEDx = []
#     ## each event should be an array of averages where each element in that array corresponds to the
#     ## respective track average at that index
#     zero_count = 0
#     for i, event_tracks in enumerate(cluster):
#         # tracks also coresponds to each event
#         #each track's average is a value
#         tracks_h1_avg = [] # array to hold all track DEDx averages in each tracks(event)
#         for track in event_tracks:
#             if len(track) == 0: #or any(x == 0 for x in track) (seems that there's no dEdx value that is 0): # no hits >> define harmonic mean = 0
#                 tracks_h1_avg.append(0)
#                 zero_count += 1
#             else:   
#                 track_avg = (len(track))/sum(x**(-1) for x in track)
#                 tracks_h1_avg.append(track_avg)
#         h1_avg_DEDx.append(tracks_h1_avg)
#     return h1_avg_DEDx
