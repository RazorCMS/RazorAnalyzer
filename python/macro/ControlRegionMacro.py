import sys
import ROOT as rt
import macro
from cutsweights import *

LUMI = 40 #in /pb
MCLUMI = 1 

def Filenames():
    return {
            "TTJets"   :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
            "WJets"    :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root",
            "SingleTop":"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/SingleTop_1pb_weighted.root",
            "Data"     :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/SingleMuonAndElectron_Run2015B-GOLDEN.root"
            }

def Fill(event, hists, weight, debug=False):
    """Fill MR and Rsq histograms"""
    varList = {"MR":event.MR, "Rsq":event.Rsq}
    for name,var in varList.iteritems(): macro.fillHist(hists, name, var, weight, debug)

def Print(histDict, c):
    """Makes plots of MR and Rsq"""
    #MC histograms
    ordering = ["TTJets", "WJets", "SingleTop"]
    for name in ordering: 
        for var in histDict[name]: macro.setHistColor(histDict[name][var], name)
    names = {name:name for name in ordering}
    MRHists = {name:histDict[name]["MR"] for name in ordering}
    RsqHists = {name:histDict[name]["Rsq"] for name in ordering}

    #data histograms
    dataName = "Data"
    dataHists = histDict[dataName]

    mrStack = macro.makeStack(MRHists, ordering, "MR")
    rsqStack = macro.makeStack(RsqHists, ordering, "Rsq")
    legend = macro.makeLegend(MRHists, names, reversed(ordering))
    legend.AddEntry(dataHists["MR"], "Data")
    macro.plot_basic(c, mc=mrStack, data=dataHists["MR"], leg=legend, xtitle="MR", printstr="MR_TTJetsSingleLepton")
    macro.plot_basic(c, mc=rsqStack, data=dataHists["Rsq"], leg=legend, xtitle="Rsq", printstr="Rsq_TTJetsSingleLepton")

if __name__ == "__main__":
    rt.gROOT.SetBatch()
    if len(sys.argv) > 1: debug = bool(sys.argv[1])
    else: debug = False
    if len(sys.argv) > 2: maxEvents = int(sys.argv[2])
    else: maxEvents = -1

    #setup
    print("Initializing")
    inputs = Filenames()
    files = {name:rt.TFile.Open(inputs[name]) for name in inputs} #get input files
    for name in inputs: assert files[name] #check input files
    hists = {name:{} for name in inputs}
    for name in inputs: 
        files[name] = rt.TFile.Open(inputs[name])
    #load trees
    trees = macro.makeTreeDict(files, "ControlSampleEvent", debug)
    #load weight histograms
    weightHists = macro.loadWeightHists(debug)

    #define histograms to fill
    for name in hists:
        hists[name]["MR"] = rt.TH1F("MR"+name, "M_{R} (GeV)", 20, 300, 4000)
        hists[name]["Rsq"] = rt.TH1F("Rsq"+name, "R^{2}", 20, 0.15, 1.5)
        for var in hists[name]: 
            hists[name][var].Sumw2()
    
    #fill histograms
    print("MC:")
    macro.loopTrees(trees, cuts_TTJetsSingleLepton, weight_mc, weightHists, Fill, {name:hists[name] for name in hists if name != "Data"}, maxEvents=maxEvents, scale=LUMI*1.0/MCLUMI, debug=debug) #MC
    print("Data:")
    macro.loopTree(trees["Data"], cuts_TTJetsSingleLepton, weight_data, {}, Fill, hists["Data"], debug=debug) #Data

    #print histograms
    c = rt.TCanvas("c", "c", 800, 600)
    Print(hists, c)

    #close files
    print("Cleaning up")
    for f in files: files[f].Close()
