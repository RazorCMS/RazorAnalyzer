import sys, os
import string
import numpy as np
import itertools
import ROOT as rt
import rootTools

##initialization

#important values
mu = 1.11 #best-fit signal strength for Higgs from Run 1
#TODO: multiply weights of higgs MC by mu

#check that the input files were specified
if len(sys.argv) < 2:
    print 'Usage: HggRazorFits.py filename1.root ...'
    exit()

#switch ROOT to batch mode
rt.gROOT.SetBatch()
#turn on fit stats
rt.gStyle.SetOptFit(0111)

#get the path of this script and load the razor library
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname) #get the full path of this script
outpath = fullpath+"/output"
libpath = fullpath+"/lib"
rt.gSystem.Load(libpath+"/libRazorRun2.so")

#suppress info messages from ROOT and RooFit
rt.RooMsgService.instance().setStreamStatus(1,False)
rt.gErrorIgnoreLevel = rt.kWarning

print('Will compute H->gg analysis signal regions and find signal yields in each box')

#load the TTrees from the input files
boxNames = [
        'HighPt',
        'Hbb',
        'Zbb',
        'HighRes',
        'LowRes'
        ]
filenames = [sys.argv[i] for i in range(1, len(sys.argv))]
files = dict((filename, rt.TFile(filename)) for filename in filenames)
#divide inputs into SM Higgs MC, signal MC, and data
higgsMCBoxes = dict((box, [file.Get(box) for filename, file in files.iteritems() if 'HToGammaGamma' in filename]) for box in boxNames)
signalMCBoxes = dict((box, [file.Get(box) for filename, file in files.iteritems() if 'SMS' in filename]) for box in boxNames)
smMCBoxes = dict((box, [file.Get(box) for filename, file in files.iteritems() if 'HToGammaGamma' not in filename and 'SMS' not in filename]) for box in boxNames)

#create RooDataSets
vars = {}
vars['mGammaGamma'] = rt.RooRealVar("mGammaGamma", "mGammaGamma", 100, 180)
vars['MR'] = rt.RooRealVar("MR", "MR", 200, 3000)
vars['Rsq'] = rt.RooRealVar("Rsq", "Rsq", 0.0, 1.0)
weight = rt.RooRealVar("weight", "weight", 0.0, 10000000.0)
args = rt.RooArgSet()
argsWeighted = rt.RooArgSet()
for key, var in vars.iteritems(): 
    args.add(var)
    argsWeighted.add(var)
argsWeighted.add(weight)

#higgs datasets
higgsMCDataSets = dict((box, rt.RooDataSet("higgsMCDataSet"+box, "Combined Higgs MC, Box "+box, argsWeighted, "weight")) for box in boxNames)
for box, trees in higgsMCBoxes.iteritems(): 
    for index, tree in enumerate(trees):
        higgsMCDataSets[box].append(rt.RooDataSet("higgsMC"+box+str(index), "higgsMC"+box+str(index), tree, argsWeighted, "", "weight"))
   
#signal datasets
signalMCDataSets = dict((box, rt.RooDataSet("signalMCDataSet"+box, "Combined Signal MC, Box "+box, argsWeighted, "weight")) for box in boxNames)
for box, trees in signalMCBoxes.iteritems():
    for index, tree in enumerate(trees):
        signalMCDataSets[box].append(rt.RooDataSet("signalMC"+box+str(index), "signalMC"+box+str(index), tree, argsWeighted, "", "weight"))

#SM datasets
smMCDataSets = dict((box, rt.RooDataSet("smMCDataSet"+box, "SM MC, Box "+box, argsWeighted, "weight")) for box in boxNames)
for box, trees in smMCBoxes.iteritems():
    for index, tree in enumerate(trees):
        smMCDataSets[box].append(rt.RooDataSet("smMC"+box+str(index), "smMC"+box+str(index), tree, argsWeighted, "", "weight"))

##step 1: use the SM Higgs MC to compute the size of the signal region
print("Determining effective Higgs peak width in each box...")
sigmaEff = dict((box, 0.0) for box in boxNames)
higgsMCSignalRegionYield = {}
higgsMCBackgroundRegionYield = {}
for box in boxNames:
    print("Finding width in "+box+" box")
    fractionHiggsEventsInSignalRegion = 0.0
    stepSize = 0.001 #step size used to widen window until 68.2% of higgs events are within
    while fractionHiggsEventsInSignalRegion < 0.682:
        sigmaEff[box] += stepSize
        fractionHiggsEventsInSignalRegion = higgsMCDataSets[box].sumEntries("mGammaGamma > "+str(125-sigmaEff[box])+" && mGammaGamma < "+str(125+sigmaEff[box]))/higgsMCDataSets[box].sumEntries()
    higgsMCSignalRegionYield[box] = higgsMCDataSets[box].sumEntries("mGammaGamma > "+str(125-2*sigmaEff[box])+" && mGammaGamma < "+str(126+2*sigmaEff[box]))
    higgsMCBackgroundRegionYield[box] = higgsMCDataSets[box].sumEntries("(mGammaGamma > 103 && mGammaGamma < 120) || (mGammaGamma > 131 && mGammaGamma < 160)")

##step 2: perform a fit in the full region 103 < mGammaGamma < 160 to extract the background yields in the signal regions

#create background PDF (sum of two exponentials)
pars = {}
pars["alpha1"] = rt.RooRealVar("alpha1", "alpha1", -0.02, -1.0, 0.0)
pars["alpha2"] = rt.RooRealVar("alpha2", "alpha2", -0.0004, -1.0, 0.0)
pars["f"] = rt.RooRealVar("f", "f", 0.005, 0.995)
pars["nEvents"] = rt.RooRealVar("N", "N", 0.0, 10000000)
expo1 = rt.RooExponential("expo1", "expo1", vars["mGammaGamma"], pars["alpha1"])
expo2 = rt.RooExponential("expo2", "expo2", vars["mGammaGamma"], pars["alpha2"])
backgroundPdf = rt.RooAddPdf("backgroundPdf", "backgroundPdf", expo1, expo2, pars["f"])
extBackgroundPdf = rt.RooExtendPdf("extBackgroundPdf", "extBackgroundPdf", backgroundPdf, pars["nEvents"])

#fit background pdf in each box and integrate to estimate number of background events in signal region
vars['mGammaGamma'].setRange('background', 103, 160)
signalRegionIntegral = {}
signalRegionIntegralErr = {}
scaleFactor = {}
scaleFactorErr = {}
backgroundCutString = "(mGammaGamma > 103 && mGammaGamma < 120) || (mGammaGamma > 131 && mGammaGamma < 160)"
c = rt.TCanvas("c", "c", 800, 600)
for box in boxNames:
    pars["alpha1"].setVal(0.3)
    pars["alpha2"].setVal(0.6)
    pars["f"].setVal(0.2)
    pars["nEvents"].setVal(smMCDataSets[box].sumEntries())
    fitResult = extBackgroundPdf.fitTo(smMCDataSets[box], rt.RooFit.Save(), rt.RooFit.Extended(rt.kTRUE), rt.RooFit.SumW2Error(rt.kTRUE))
    fitResult.Print("v")
    vars["mGammaGamma"].setRange("signal"+box, 125-2*sigmaEff[box], 126+2*sigmaEff[box])
    signalRegionIntegralObj = extBackgroundPdf.createIntegral(rt.RooArgSet(vars['mGammaGamma']), rt.RooFit.Range('signal'+box))
    signalRegionIntegral[box] = signalRegionIntegralObj.getVal()
    signalRegionIntegralErr[box] = 0 #TODO: find out how to get this error

    #plot the fit result
    frame = vars["mGammaGamma"].frame()
    smMCDataSets[box].plotOn(frame)
    extBackgroundPdf.plotOn(frame)
    frame.SetTitle("Background fit in "+box+" box")
    frame.Draw()
    c.Print(outpath+"/HggBackgroundFit"+box+".pdf")

    #compute background prediction scale factor: integral of background pdf in sig region / actual num in sideband
    #background region: 103 < mGammaGamma < 120 or 131 < mGammaGamma < 160
    numInSideband = smMCDataSets[box].sumEntries(backgroundCutString)
    if numInSideband > 0: 
        scaleFactor[box] = signalRegionIntegral[box]/numInSideband
        scaleFactorErr[box] = signalRegionIntegralErr[box]/numInSideband #TODO: include error on sideband yield
    else: 
        print("Number of events in sideband is zero!  Setting scale factor for box "+box+" to 0.0.")
        scaleFactor[box] = 0.0
        scaleFactorErr[box] = 0.0

##compute signal region yields from signal MC
signalMCSignalRegionYield = {}
signalMCBackgroundRegionYield = {}
for box in boxNames:
    signalMCSignalRegionYield[box] = signalMCDataSets[box].sumEntries("mGammaGamma > "+str(125-2*sigmaEff[box])+" && mGammaGamma < "+str(126+2*sigmaEff[box]))
    signalMCBackgroundRegionYield[box] = signalMCDataSets[box].sumEntries(backgroundCutString)

##step 3: plot the distribution of sideband events in the MR-Rsq plane
combinedMCDataSets = {}
MRMaxValue = 3000
print("\nComputing signal regions...\n")
for box in boxNames:
    backgroundHist = smMCDataSets[box].createHistogram(vars["MR"], vars["Rsq"], 300, 500, backgroundCutString)
    backgroundHist.GetXaxis().SetRangeUser(200, 1000)
    backgroundHist.GetYaxis().SetRangeUser(0.0, 0.2)
    backgroundHist.GetZaxis().SetRangeUser(0.0, 200)
    backgroundHist.SetTitle("Nonresonant background distribution in "+box+" box")
    backgroundHist.SetStats(0)
    backgroundHist.Draw("colz")
    c.Print(outpath+"/HggRazorBackgroundDistribution"+box+".pdf")
    hggHist = higgsMCDataSets[box].createHistogram(vars["MR"], vars["Rsq"])
    hggHist.GetXaxis().SetRangeUser(200, 3000)
    hggHist.GetYaxis().SetRangeUser(0.0, 1)
    hggHist.GetZaxis().SetRangeUser(0.0, 1.0)
    hggHist.SetTitle("Higgs background distribution in "+box+" box")
    hggHist.SetStats(0)
    hggHist.Draw("colz")
    c.Print(outpath+"/HggRazorHiggsDistribution"+box+".pdf")
    combinedMCDataSets[box] = smMCDataSets[box].Clone("combinedMCDataSets"+box)
    combinedMCDataSets[box].append(higgsMCDataSets[box])
    #create signal regions in the MR-Rsq plane
    RsqCut = 0
    MRCut = MRMaxValue
    MREdges = [MRMaxValue]
    done = False
    while not done:
        #compute number of expected events in the region [MRCut, minMREdge]x[RsqCut, 1]
        Integral = combinedMCDataSets[box].sumEntries("Rsq < 1 && Rsq > "+str(RsqCut)+" && MR > "+str(MRCut)+" && MR < "+str(min(MREdges)))
        #if there is at least one event, make a new MR edge and increment the Rsq cut
        if Integral >= 1:
            MREdges.append(MRCut)
            RsqCut += 0.05
        #if the MR cut is less than 150, add it and finish
        if MRCut <= 150: 
            done = True
            MREdges.append(MRCut)
        #decrease the MR cut for the next round
        if MRCut < 500: MRCut = MRCut - 50
        elif MRCut < 1000: MRCut = MRCut - 100
        else: MRCut = MRCut - 200
    MREdges.remove(MRMaxValue) #for convenience in creating the bins 
    #create the bins defined by the above algorithm
    signalBinCorners = []
    MRHigh = MRMaxValue
    for index, edge in enumerate(MREdges):
        MRLow = edge
        RsqHigh = 0.0
        for i in range(index):
            RsqLow = RsqHigh
            RsqHigh += 0.05
            signalBinCorners.append([[MRLow, MRHigh], [RsqLow, RsqHigh]])
            #make a bin
        #make the last bin for this MR value
        signalBinCorners.append([[MRLow, MRHigh], [RsqHigh, 1.0]])
        #set max MR for next bin
        MRHigh = MRLow
    #print out yields in each signal region
    print("Box "+box)
    print("     Yields in box:")
    print("         SM MC:     "+str(smMCDataSets[box].sumEntries()))
    print("         Higgs MC:  "+str(higgsMCDataSets[box].sumEntries()))
    print("         Signal MC: "+str(signalMCDataSets[box].sumEntries()))

    print("     SigmaEff = "+str(sigmaEff[box]))
    print("     Integral of background pdf from "+str(125-2*sigmaEff[box])+" to "+str(126+2*sigmaEff[box])+" = "+str(signalRegionIntegral[box]))
    print("     Scale factor = "+str(scaleFactor[box])+" +/- "+str(scaleFactorErr[box]))
    print("     Events in signal region: Higgs MC "+str(higgsMCSignalRegionYield[box])+", Signal MC "+str(signalMCSignalRegionYield[box]))
    print("     Events in background region: Higgs MC "+str(higgsMCBackgroundRegionYield[box])+", Signal MC "+str(signalMCBackgroundRegionYield[box]))
    print("     MR bin edges: "+str(MREdges)+" and 3000")
    print("     Yields by signal region: ")
    signalRegionYieldsMC = []
    for [[MRMin, MRMax], [RsqMin, RsqMax]] in signalBinCorners:
        cutString = "MR > "+str(MRMin)+" && MR < "+str(MRMax)+" && Rsq > "+str(RsqMin)+" && Rsq < "+str(RsqMax)
        nSig = combinedMCDataSets[box].sumEntries(cutString)
        signalRegionYieldsMC.append(nSig)
        print("         "+str([[MRMin, MRMax], [RsqMin, RsqMax]])+": "+str(nSig))
