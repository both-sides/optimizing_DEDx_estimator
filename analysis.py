import ROOT as rt
import numpy as np


f = rt.TFile.Open("data/HSCPgluino_M-1800_fromAOD.root")
dir = f.GetDirectory("HSCPFullAODAnalyzer")
tree = f.Get("HSCPFullAODAnalyzer/Events")

filename = "data/HSCPgluino_M-1800_fromAOD.root"
treename = "HSCPFullAODAnalyzer/" + dir.Get("Events").GetName()

df = rt.RDataFrame(treename, filename)

branches = [b.GetName() for b in tree.GetListOfBranches()]

with open("branches.txt", "w") as ff:
    for branch in branches:
        ff.write(branch + "\n")

f.Close()

