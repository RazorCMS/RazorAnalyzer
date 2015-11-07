### fit MC cocktail to validate the fit function

import sys,os
import argparse
import ROOT as rt

#local imports
from framework import Config
import macro.macro as macro
from macro.razorAnalysis import *
from macro.razorWeights import *
from macro.razorMacros import *

LUMI = 32000 #in /pb
MCLUMI = 1 

SAMPLES = ["TTV", "VV", "DYJetsLow", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]

DIR_MC= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC"
PREFIX = "RazorInclusive"
FILENAMES = {
        "TTJets"    : DIR_MC+"/"+PREFIX+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted_razorskim.root",
        "WJets"     : DIR_MC+"/"+PREFIX+"_WJetsToLNu_HTBinned_1pb_weighted_razorskim.root",
        "SingleTop" : DIR_MC+"/"+PREFIX+"_SingleTop_1pb_weighted_razorskim.root",
        "DYJetsLow" : DIR_MC+"/"+PREFIX+"_DYJetsToLL_M-5to50_HTBinned_1pb_weighted.root",
        "VV" : DIR_MC+"/"+PREFIX+"_VV_1pb_weighted_razorskim.root",
        "TTV" : DIR_MC+"/"+PREFIX+"_TTV_1pb_weighted_razorskim.root",
        "DYJets"     : DIR_MC+"/"+PREFIX+"_DYJetsToLL_M-50_HTBinned_1pb_weighted_razorskim.root",
        "ZInv"     : DIR_MC+"/"+PREFIX+"_ZJetsToNuNu_HTBinned_1pb_weighted_razorskim.root",
        }

config = "config/run2_sideband.config"
FIT_DIR = "FitPlots"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/toys_Bayes_MultiJet.root",
        "MuMultiJet":FIT_DIR+"/toys_Bayes_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/toys_Bayes_EleMultiJet.root",
        }

cfg = Config.Config(config)
binsMRHad = cfg.getBinning("MultiJet")[0]
binsRsqHad = cfg.getBinning("MultiJet")[1]
hadronicBinning = { "MR":binsMRHad, "Rsq":binsRsqHad }
binsMRLep = cfg.getBinning("MuMultiJet")[0]
binsRsqLep = cfg.getBinning("MuMultiJet")[1]
leptonicBinning = { "MR":binsMRLep, "Rsq":binsRsqLep }
binning = { "MultiJet":hadronicBinning, "MuMultiJet":leptonicBinning, "EleMultiJet":leptonicBinning}

weightOpts = []
shapeErrors = []
miscErrors = []

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
    weightHists = loadWeightHists(weightfilenames_DEFAULT, weighthistnames_DEFAULT, debugLevel)
    sfHists = {}

    #estimate yields in leptonic signal region
    for lepType in ["", "Mu", "Ele"]:
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
                    runFitAndToys(FIT_DIR, boxName, LUMI, PREFIX+'_'+DATA_NAME, DIR_DATA, config=config, sideband=True)
                    #check
                    if not os.path.isfile(TOYS_FILES[boxName]):
                        print "Error creating fit file",TOYS_FILES[boxName]
                        sys.exit()
                makeControlSampleHists(extboxName, 
                        filenames=FILENAMES, samples=SAMPLES, 
                        cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                        bins=binning[boxName], lumiMC=MCLUMI, lumiData=LUMI, 
                        weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                        weightOpts=weightOpts, shapeErrors=shapeErrors, miscErrors=miscErrors,
                        fitToyFiles=TOYS_FILES, boxName=boxName, 
                        btags=nBtags, debugLevel=debugLevel)

