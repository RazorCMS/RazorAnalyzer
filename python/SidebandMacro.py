import sys,os
import argparse
import ROOT as rt

#local imports
import macro.macro as macro
from macro.razorAnalysis import *
from macro.razorWeights import *
from macro.razorMacros import *

LUMI = 1264 #in /pb
MCLUMI = 1 

SAMPLES_SIGNAL = ["TTV", "VV", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]

DIR_MC= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC_FullRazorInclusive/combined"
DIR_DATA= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data"
PREFIX_SIGNAL = "RazorInclusive"
DATA_NAME='SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_razorskim'
FILENAMES_SIGNAL = {
        "TTJets"    : DIR_MC+"/"+PREFIX_SIGNAL+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root",
        "WJets"     : DIR_MC+"/"+PREFIX_SIGNAL+"_WJets_1pb_weighted.root",
        "SingleTop" : DIR_MC+"/"+PREFIX_SIGNAL+"_SingleTop_1pb_weighted.root",
        "VV" : DIR_MC+"/"+PREFIX_SIGNAL+"_VV_1pb_weighted.root",
        "TTV" : DIR_MC+"/"+PREFIX_SIGNAL+"_TTV_1pb_weighted.root",
        "DYJets"     : DIR_MC+"/"+PREFIX_SIGNAL+"_DYJetsToLL_1pb_weighted.root",
        "ZInv"     : DIR_MC+"/"+PREFIX_SIGNAL+"_ZJetsToNuNu_1pb_weighted.root",
        "Data"      : DIR_DATA+'/'+PREFIX_SIGNAL+'_'+DATA_NAME+'.root',
        }

config = "config/backup.config"
FIT_DIR = "FitResultsSideband"
TOYS_FILES = {
        "MuMultiJet":FIT_DIR+"/toys_Bayes_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/toys_Bayes_EleMultiJet.root",
        }

WEIGHTDIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors"
LEPTONWEIGHTDIR = "LeptonEfficiencies/20151013_PR_2015D_Golden_1264"
weightfilenames = {
        "muon": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D_Golden.root",
        "ele": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root",
        "muontrig": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015D_Golden.root",
        "eletrig": WEIGHTDIR+"/"+LEPTONWEIGHTDIR+"/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015D_Golden.root",
        "pileup": WEIGHTDIR+"/PileupWeights/NVtxReweight_ZToMuMu_2015D_1264ipb.root",
        }
weighthistnames = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "muontrig": "ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight",
        "eletrig": "ScaleFactor_EleTriggerEleCombinedEffDenominatorTight",
        "pileup": "NVtxReweight",
        }
weightOpts = ["doNVtxWeights"]
shapeErrors = ["muoneff", "eleeff", "jes"]
miscErrors = ["mt"]

oldBlindBins = [(x,y) for x in range(2,len(leptonicSidebandBins["MR"])+1) for y in range(2,len(leptonicSidebandBins["Rsq"])+1)]

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
    weightHists = loadWeightHists(weightfilenames, weighthistnames, debugLevel)
    sfHists = {}

    #get scale factor histograms
    sfHists = loadScaleFactorHists(processNames=SAMPLES_SIGNAL, debugLevel=debugLevel)

    #estimate yields in signal region
    for lepType in ["Mu", "Ele"]:
        for jets in ["MultiJet"]:
            boxName = lepType+jets
            #btaglist = [0]
            btaglist = [0,1,2,3]
            for btags in btaglist:
                print "\n---",boxName,"Box,",btags,"B-tags ---"
                #get correct cuts string
                thisBoxCuts = razorCuts[boxName]
                if btags < len(btaglist)-1:
                    thisBoxCuts += " && nBTaggedJets == "+str(btags)
                else:
                    thisBoxCuts += " && nBTaggedJets >= "+str(btags)

                if len(btaglist) > 1:
                    extboxName = boxName+str(btags)+"BTag"
                    nBtags = btags
                else:
                    extboxName = boxName
                    nBtags = -1
                #check fit file and create if necessary
                if not os.path.isfile(TOYS_FILES[boxName]):
                    print "Fit file",TOYS_FILES[boxName],"not found, trying to recreate it"
                    runFitAndToys(FIT_DIR, boxName, LUMI, PREFIX_SIGNAL+'_'+DATA_NAME, DIR_DATA, config=config, sideband=True)
                    #check
                    if not os.path.isfile(TOYS_FILES[boxName]):
                        print "Error creating fit file",TOYS_FILES[boxName]
                        sys.exit()
                makeControlSampleHists(extboxName, 
                        filenames=FILENAMES_SIGNAL, samples=SAMPLES_SIGNAL, 
                        cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                        bins=leptonicSignalRegionBins, lumiMC=MCLUMI, lumiData=LUMI, 
                        weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                        weightOpts=weightOpts, shapeErrors=shapeErrors, miscErrors=miscErrors,
                        fitToyFiles=TOYS_FILES, boxName=boxName, blindBins=oldBlindBins,
                        btags=nBtags, debugLevel=debugLevel)
