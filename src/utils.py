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

