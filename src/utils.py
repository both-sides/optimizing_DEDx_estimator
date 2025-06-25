import math
import ROOT as rt
from typing import List
from datetime import date


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
    13: rt.kGray,  # neutral mid‐tone
    # two extra custom colors:
    14: rt.TColor.GetColor("#8B4513"),  # Brown
    15: rt.TColor.GetColor("#00CED1"),  # DarkTurquoise
}


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
    
    
    #Two‐panel ratio plot: top = overlay, bottom = ratio.
    def draw_ratio(self,
                   h_num: rt.TH1,
                   h_den: rt.TH1,
                   canvas_name: str = "ratio",
                   title: str = "",
                   xlabel: str = "",
                   ylabel: str = "",
                   ratio_ylabel: str = "Ratio",
                   ratio_range: tuple[float,float] = (-5, 5),
                   logy: bool = False):
        
        # close any existing canvases so ROOT will make a truly new one
        for canv in rt.gROOT.GetListOfCanvases():
            canv.Close()
    
        # 1) make your canvas
        c = rt.TCanvas("c_ratio","Ratio Plot", 800, 800)

        #create two pads: top (0.3–1) and bottom (0–0.3)
        pad1 = rt.TPad(f"{canvas_name}_pad1", "Top pad", 0, 0.3, 1, 1)
        pad2 = rt.TPad(f"{canvas_name}_pad2","Bottom pad", 0, 0,   1, 0.3)
        #pad1.SetBottomMargin(0)   # no x-axis tick labels
        pad2.SetTopMargin(0.02)      # no title overlap
        pad2.SetBottomMargin(0.4)    # room for x-axis label
        pad1.Draw() 
        pad2.Draw()

        # 2) top pad: overlay
        pad1.cd()
        if logy: 
            pad1.SetLogy(1)
        h_den.SetTitle(title)
        h_den.SetLineColor(rt.kBlack)
        h_num.SetLineColor(rt.kRed)
        
        # auto-scale so they both fit
        maxval = max(h_den.GetMaximum(), h_num.GetMaximum())
        h_den.SetMaximum(maxval*1.2)
        h_den.SetTitle(title)
        h_den.GetYaxis().SetTitle("Number of tracks")
        h_den.Draw()
        h_num.Draw("same")
        pad1.Update()


        #build the ratio histogram in pad2
        pad2.cd()
        ratio = h_num.Clone("ratio")       # clone numerator
        ratio.Divide(h_den)                # divide by denominator
        ratio.SetMarkerStyle(20)           # draw as points
        ratio.SetTitle("")                 # no global title
        ratio.GetYaxis().SetTitle("Num/Den")
        ratio.GetXaxis().SetTitle("DEDx_IhStrip (MeV/cm)")
        # tweak axis label sizes so they’re big in the small pad:
        ratio.GetYaxis().SetTitleSize(0.05)
        ratio.GetYaxis().SetTitleOffset(0.4)
        ratio.GetYaxis().SetLabelSize(0.05)
        ratio.GetXaxis().SetTitleSize(0.08)
        ratio.GetXaxis().SetLabelSize(0.08)

        # set a sensible ratio range, e.g. 0.5–1.5
        ratio.SetMinimum(-5)
        ratio.SetMaximum(5)
        ratio.Draw("EP")                   # E: error bars, P: points
        pad2.Update()


        #finally draw the overall canvas
        c.cd()
        c.Modified()
        c.Draw()
        return c
                   
        
        
                              
        
                   
    
    
                 

