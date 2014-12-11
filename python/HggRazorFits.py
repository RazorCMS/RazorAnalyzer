import sys, os
import string
import numpy as np
import itertools
import ROOT as rt
import rootTools

##initialization

#important values
integratedLumi = 5000 #in /pb
mu = 1.11 #best-fit signal strength for Higgs from Run 1

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

#suppress info messages from RooFit
#rt.RooMsgService.instance().setStreamStatus(1,False)

print('Will compute H->gg analysis signal regions and compute signal yields in each box')

#load the TTrees from the input files
boxNames = [
        'HighPt',
        'Hbb',
        'Zbb',
        'HighRes',
        'LowRes'
        ]
filenames = [sys.argv[i] for i in range(1, len(sys.argv))]
files = [rt.TFile(filename) for filename in filenames]
#divide inputs into SM Higgs MC, signal MC, and data
higgsMCBoxes = dict((box, [file.Get(box) for file in files]) for box in boxNames if 'HToGG' in box)
signalMCBoxes = dict((box, [file.Get(box) for file in files]) for box in boxNames if 'SMS' in box)
smMCBoxes = dict((box, [file.Get(box) for file in files]) for box in boxNames if 'SMS' not in box and 'HToGG' not in box)

#create RooDataSets
vars = {}
vars['mGammaGamma'] = rt.RooRealVar("mGammaGamma", "mGammaGamma", 100, 180)
vars['MR'] = rt.RooRealVar("MR", "MR", 200, 4000)
vars['Rsq'] = rt.RooRealVar("Rsq", "Rsq", 0.0, 2.0)
weight = rt.RooRealVar("weight", "weight", 0.0, 10000000.0)
args = rt.RooArgSet()
argsWeighted = rt.RooArgSet()
for var in vars: 
    args.add(var)
    argsWeighted.add(var)
argsWeighted.add(weight)

#higgs datasets
higgsMCDataSetsOld = dict((box, rt.RooDataSet("higgsMCDataSetOld"+box, "Combined Higgs MC, Box "+box, argsWeighted, "weight")) for box in boxNames if 'HToGG' in box)
higgsMCDataSets = dict((box, rt.RooDataSet("higgsMCDataSet"+box, "Combined Higgs MC, Box "+box, argsWeighted, "weight")) for box in boxNames if 'HToGG' in box) #with weights multiplied by integrated lumi * mu
for (box, trees) in higgsMCBoxes: 
    for index, tree in enumerate(trees):
        higgsMCDataSetsOld[box].append(rt.RooDataSet("higgsMC"+box+str(index), "higgsMC"+box+str(index), tree, argsWeighted, "", "weight"))
    #rescale weights
    for i in range(0, higgsMCDataSetsOld[box].numEntries()):
        theseArgs = higgsMCDataSetOld[box].get(i)
        theseArgs.setRealValue("weight", theseArgs.getRealValue("weight")*integratedLumi*mu)
        higgsMCDataSets[box].add(theseArgs)
    
#signal datasets
signalMCDataSets = dict((box, rt.RooDataSet("signalMCDataSet"+box, "Combined Signal MC, Box "+box, argsWeighted, "weight")) for box in boxNames if 'SMS' in box)
signalMCDataSetsOld = dict((box, rt.RooDataSet("signalMCDataSetOld"+box, "Combined Signal MC, Box "+box, argsWeighted, "weight")) for box in boxNames if 'SMS' in box)
for (box, trees) in signalMCBoxes:
    for index, tree in enumerate(trees):
        signalMCDataSetsOld[box].append(rt.RooDataSet("signalMC"+box+str(index), "signalMC"+box+str(index), tree, argsWeighted, "", "weight"))
    #rescale weights
    for i in range(0, signalMCDataSetsOld[box].numEntries()):
        theseArgs = signalMCDataSetOld[box].get(i)
        theseArgs.setRealValue("weight", theseArgs.getRealValue("weight")*integratedLumi)
        signalMCDataSets[box].add(theseArgs)

#SM datasets
smMCDataSetsOld = dict((box, rt.RooDataSet("smMCDataSetOld"+box, "Combined Data, Box "+box, argsWeighted, "weight")) for box in boxNames if 'HToGG' not in box and 'SMS' not in box)
smMCDataSets = dict((box, rt.RooDataSet("smMCDataSet"+box, "Combined Data, Box "+box, argsWeighted, "weight")) for box in boxNames if 'HToGG' not in box and 'SMS' not in box)
for (box, trees) in smMCBoxes:
    for index, tree in enumerate(trees):
        smMCDataSetsOld[box].append(rt.RooDataSet("smMC"+box+str(index), "smMC"+box+str(index), tree, argsWeighted, "", "weight"))
    #rescale weights
    for i in range(0, smMCDataSetsOld[box].numEntries()):
        theseArgs = smMCDataSetOld[box].get(i)
        theseArgs.setRealValue("weight", theseArgs.getRealValue("weight")*integratedLumi)
        smMCDataSets[box].add(theseArgs)

##step 1: use the SM Higgs MC to compute the size of the signal region
print("Determining effective Higgs peak width in each box...")
sigmaEff = dict((box, 0.0) for box in boxNames)
higgsMCSignalRegionYield = {}
higgsMCBackgroundRegionYield = {}
for box in boxNames:
    fractionHiggsEventsInSignalRegion = 0.0
    stepSize = 0.001 #step size used to widen window until 68.2% of higgs events are within
    while fractionHiggsEventsInSignalRegion < 0.682:
        sigmaEff[box] += stepSize
        fractionHiggsEventsInSignalRegion = higgsMCDataSets[box].sumEntries("mGammaGamma > "+str(125-sigmaEff[box])+" && mGammaGamma < "+str(125+sigmaEff[box]))/higgsMCDataSets[box].sumEntries()
    print("sigmaEff in box "+box+" = "+str(sigmaEff[box]))
    higgsMCSignalRegionYield[box] = higgsMCDataSets[box].sumEntries("mGammaGamma > "str(125-2*sigmaEff[box])+" && mGammaGamma < "+str(126+2*sigmaEff[box]))
    higgsMCBackgroundRegionYield[box] = higgsMCDataSets[box].sumEntries("(mGammaGamma > 103 && mGammaGamma < 120) || (mGammaGamma > 131 && mGammaGamma < 160)")

##step 2: perform a fit in the full region 103 < mGammaGamma < 160 to extract the background yields in the signal regions

#create background PDF (sum of two exponentials)
pars = {}
pars["alpha1"] = rt.RooRealVar("alpha1", "alpha1", -1.0, 0.0)
pars["alpha2"] = rt.RooRealVar("alpha2", "alpha2", -1.0, 0.0)
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
    nEvents.setVal(smMCDataSets[box].sumEntries())
    fitResult = extBackgroundPdf.fitTo(smMCDataSets[box], rt.RooFit.Extended(rt.kTRUE), rt.RooFit.SumW2Errors(rt.kTRUE))
    fitResult.Print("v")
    vars["mGammaGamma"].setRange("signal"+box, 125-2*sigmaEff[box], 126+2*sigmaEff[box])
    signalRegionIntegralObj = extBackgroundPdf.createIntegral(vars['mGammaGamma'], rt.RooFit.Range('signal'+box))
    signalRegionIntegral[box] = signalRegionIntegralObj.getVal()
    signalRegionIntegralErr[box] = signalRegionIntegralObj.getErr()
    print("Integral of background pdf in box "+box+" from "+str(125-2*sigmaEff[box])+" to "+str(126+2*sigmaEff[box])+" = "+str(signalRegionIntegral[box]))

    #plot the fit result
    frame = vars["mGammaGamma"].frame()
    smMCDataSets[box].plotOn(frame)
    extBackgroundPdf.plotOn(frame)
    frame.Draw()
    c.Print(outpath+"/HggBackgroundFit"+box+".pdf")

    #compute background prediction scale factor: integral of background pdf in sig region / actual num in sideband
    #background region: 103 < mGammaGamma < 120 or 131 < mGammaGamma < 160
    numInSideband = smMCDataSets[box].sumEntries(backgroundCutString)
    scaleFactor[box] = signalRegionIntegral[box]/numInSideband
    scaleFactorErr[box] = signalRegionIntegralErr[box]/numInSideband #TODO: include error on sideband yield
    print("Scale factor in box "+box+" = "+str(scaleFactor[box])+" +/- "+str(scaleFactorErr[box]))

##print out signal region yields from signal MC
signalMCSignalRegionYield = {}
signalMCBackgroundRegionYield = {}
for box in boxNames:
    signalMCSignalRegionYield[box] = signalMCDataSets[box].sumEntries("mGammaGamma > "str(125-2*sigmaEff[box])+" && mGammaGamma < "+str(126+2*sigmaEff[box]))
    signalMCBackgroundRegionYield[box] = signalMCDataSets[box].sumEntries(backgroundCutString)
    print("Box "+box)
    print("Events in signal region: Higgs MC "+str(higgsMCSignalRegionYield[box])+", Signal MC "+str(signalMCSignalRegionYield[box]))
    print("Events in background region: Higgs MC "+str(higgsMCBackgroundRegionYield[box])+", Signal MC "+str(signalMCBackgroundRegionYield[box]))

##step 3: plot the distribution of sideband events in the MR-Rsq plane
combinedMCDataSets = {}
MRMaxValue = 3000
for box in boxNames:
    backgroundHist = smMCDataSet[box].createHistogram(vars["MR"], vars["Rsq"], backgroundCutString)
    backgroundHist.Draw("colz")
    c.Print(outpath+"/HggRazorBackgroundDistribution"+box+".pdf")
    hggHist = higgsMCDataSet[box].createHistogram(vars["MR"], vars["Rsq"])
    hggHist.Draw("colz")
    c.Print(outpath+"/HggRazorHiggsDistribution"+box+".pdf")
    combinedMCDataSets[box] = smMCDataSets[box]
    combinedMCDataSets[box].append(higgsMCDataSets[box])
    #create signal regions in the MR-Rsq plane
    RsqCut = 0
    MRCut = MRMaxValue
    MREdges = [MRMaxValue]
    done = False
    while(!done):
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
        else if MRCut < 1000: MRCut = MRCut - 100
        else: MRCut = MRCut - 200
    MREdges.remove(MRMaxValue) #for convenience in creating the bins 
    print("MR bin edges: "+str(MREdges)+" and 3000")
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
    print("Signal region bins created using the algorithm: "+str(signalBinCorners))
    #print out yields in each signal region
    print("Yields by signal region: ")
    signalRegionYieldsMC = []
    for [[MRMin, MRMax], [RsqMin, RsqMax]] in signalBinCorners:
        cutString = "MR > "+str(MRMin)+" && MR < "+str(MRMax)+" && Rsq > "+str(RsqMin)+" && Rsq < "+str(RsqMax)
        nSig = combinedMCDataSets[box].sumEntries(cutString)
        signalRegionYieldsMC.append(nSig)
        print(str([[MRMin, MRMax], [RsqMin, RsqMax]])+": "+str(nSig))
