import sys
import ROOT as rt
import macro
from cutsweights import *

def Fill(event, hists, weight, debug=False):
    """Fill MR and Rsq histograms"""
    varList = {}
    varList["MR"] = event.MR
    varList["Rsq"] = event.Rsq
    for name,var in varList.iteritems(): macro.fillHist(hists, name, var, weight, debug)

def Print(histDict):
    """Makes a plot of MR"""
    ordering = ["SingleTop"]

    names = {}
    MRHists = {}
    for name in histDict: names[name] = name
    for name in histDict: MRHists[name] = histDict[name]["MR"]

    stack = macro.makeStack(MRHists, ordering, "MR")
    legend = macro.makeLegend(MRHists, names, reversed(ordering))
    macro.plot_basic(mc=stack, leg=legend, xtitle="MR", printstr="MR_TTJetsSingleLepton")

if __name__ == "__main__":
    rt.gROOT.SetBatch()
    if len(sys.argv) > 1: debug = bool(sys.argv[1])
    else: debug = False
    maxEvents = 10000

    #setup
    inputs = {"SingleTop":"SingleTop_1pb_weighted.root"}
    #load files
    files = {}
    for name in inputs: 
        files[name] = rt.TFile(inputs[name])
        assert files[name]
    #load trees
    trees = macro.makeTreeDict(files, "ControlSampleEvent", debug)
    #load weight histograms
    weightHists = macro.loadWeightHists(debug)

    #define histograms to fill
    hists = {"SingleTop": {}}
    for name in hists:
        hists[name]["MR"] = rt.TH1F("MR"+name, "M_{R} (GeV)", 20, 300, 4000)
        for var in hists[name]: 
            hists[name][var].Sumw2()
            macro.setHistColor(hists[name][var], name)
    
    macro.loopTrees(trees, cuts_TTJetsSingleLepton, weights_standardRun2, weightHists, Fill, hists, maxEvents=maxEvents, debug=debug)
    Print(hists)

#TODO: add other TTJets samples
#TODO: expand to other CRs
#TODO: add fit
