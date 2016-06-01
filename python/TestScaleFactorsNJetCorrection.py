import sys, os
import argparse
import ROOT as rt

#local imports
from framework import Config
from macro import macro
from macro.razorAnalysis import OneLeptonForNJetsCutsData, OneLeptonForNJetsCutsMC, OneLeptonInvForNJetsCutsData, OneLeptonInvForNJetsCutsMC 
from macro.razorWeights import loadScaleFactorHists
from macro.razorMacros import makeControlSampleHists, appendScaleFactors
from ComputeScaleFactorsMacro import xbinsWJETS1L, colsWJETS1L, xbinsTTJETS1L, colsTTJETS1L, xbinsWJETS1LINV, colsWJETS1LINV
from SidebandMacro import LUMI as LUMI_DATA
from DoScaleFactorsNJetCorrection import MCLUMI, SAMPLES_TTJ1L, SAMPLES_WJ1L_INV, FILENAMES_1L, FILENAMES_1L_INV, ControlRegionBinning, ControlRegionBinning_Inv, ScaleFactorVars_WJETS1L_INV

weightOpts = []

printdir="ScaleFactorsNJetPlots"

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug

    #initialize
    weightHists = {}

    #make output directory
    os.system('mkdir -p '+printdir)

    sfVars = ("MR","Rsq")
    sfVars_NoW = ("MR_NoW", "Rsq_NoW")
    auxSFs = {"NJets":("NJets40","1")} 

    plotOpts = { "comment":False }

    sfHists = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root", processNames=SAMPLES_TTJ1L, debugLevel=debugLevel)
    sfNJetsFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors.root")
    sfHists['NJets'] = sfNJetsFile.Get("NJetsCorrectionScaleFactors")

    #define tests
    tests = ["OneLeptonDiJetClosureTest"]
    cutsData = [OneLeptonForNJetsCutsData+' && NJets80 > 1 && NJets40 < 4']
    cutsMC = [OneLeptonForNJetsCutsMC+' && NJets80 > 1 && NJets40 < 4']
    btags = [-1]
    tests.append("OneLeptonJetClosureTest")
    cutsData.append(OneLeptonForNJetsCutsData+' && NJets40 < 4')
    cutsMC.append(OneLeptonForNJetsCutsMC+' && NJets40 < 4')
    btags.append(-1)
    for nb in range(4):
        tests.append("OneLeptonDiJet"+str(nb)+"BJetClosureTest")
        cutsData.append(OneLeptonForNJetsCutsData+' && NJets80 > 1 && NJets40 < 4 && NBJetsMedium =='+str(nb))
        cutsMC.append(OneLeptonForNJetsCutsMC+' && NJets80 > 1 && NJets40 < 4 && NBJetsMedium =='+str(nb))
        btags.append(nb)
        tests.append("OneLeptonJet"+str(nb)+"BJetClosureTest")
        cutsData.append(OneLeptonForNJetsCutsData+' && NJets40 < 4 && NBJetsMedium =='+str(nb))
        cutsMC.append(OneLeptonForNJetsCutsMC+' && NJets40 < 4 && NBJetsMedium =='+str(nb))
        btags.append(nb)
    #tests.append("OneLeptonMultiJetClosureTest")
    #cutsData.append(OneLeptonForNJetsCutsData+' && NJets40 > 3')
    #cutsMC.append(OneLeptonForNJetsCutsMC+' && NJets40 > 3')
    #btags.append(-1)
    #for nb in range(4):
    #    tests.append("OneLeptonMultiJet"+str(nb)+"BJetClosureTest")
    #    cutsData.append(OneLeptonForNJetsCutsData+' && NJets40 > 3 && NBJetsMedium =='+str(nb))
    #    cutsMC.append(OneLeptonForNJetsCutsMC+' && NJets40 > 3 && NBJetsMedium =='+str(nb))
    #    btags.append(nb)

    ##########################################################
    #1L W+Jets & TTJets control sample (2 jets)
    ##########################################################    
    hists = {}
    for i,t in enumerate(tests):
        hists[t] = makeControlSampleHists(t, 
                      filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                      cutsMC=cutsMC[i], cutsData=cutsData[i], 
                      bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, auxSFs=auxSFs,
                      weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, 
                      printdir=printdir, btags=btags[i], plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
                      plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
