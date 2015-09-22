import sys
import argparse
import ROOT as rt

#local imports
import macro
from razorAnalysis import *

LUMI = 40 #in /pb
MCLUMI = 1 

SAMPLES = ["TTJets", "WJets", "SingleTop"]
FILENAMES = {
            "TTJets"   :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
            "WJets"    :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root",
            "SingleTop":"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/SingleTop_1pb_weighted.root",
            "Data"     :"root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/V1p16/OneLeptonReduced_new2/SingleMuonAndElectron_Run2015B-GOLDEN.root"
            }

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="display detailed output messages",
                                action="store_true")
    args = parser.parse_args()
    debug = args.debug

    #setup files and trees
    inputs = FILENAMES
    files = {name:rt.TFile.Open(inputs[name]) for name in inputs} #get input files
    for name in inputs: assert files[name] #check input files
    trees = macro.makeTreeDict(files, "ControlSampleEvent", debug)

    #load weight histograms
    weightHists = loadWeightHists(debug)

    #define histograms to fill
    hists = {name:{} for name in inputs}
    for name in inputs:
        hists[name]["MR"] = rt.TH1F("MR"+name, "M_{R} (GeV)", 20, 300, 4000)
        hists[name]["Rsq"] = rt.TH1F("Rsq"+name, "R^{2}", 20, 0.15, 1.5)
        for var in hists[name]: 
            hists[name][var].Sumw2()
    listOfVars = hists.itervalues().next().keys() #list of the variable names
    
    #fill histograms
    print("MC:")
    macro.loopTrees(trees, weightF=weight_mc, cuts=ttjetsSingleLeptonCutsMC, varList=listOfVars, hists={name:hists[name] for name in SAMPLES}, weightHists=weightHists, scale=LUMI*1.0/MCLUMI, debug=debug) 

    print("Data:")
    macro.loopTree(trees["Data"], weightF=weight_data, cuts=ttjetsSingleLeptonCutsData, varList=hists["Data"].keys(), hists=hists["Data"], weightHists=weightHists, debug=debug) 

    #print histograms
    c = rt.TCanvas("c", "c", 800, 600)
    macro.basicPrint(hists, mcNames=SAMPLES, varList=listOfVars, c=c, printName="TTJetsSingleLepton")

    #close files
    for f in files: files[f].Close()
