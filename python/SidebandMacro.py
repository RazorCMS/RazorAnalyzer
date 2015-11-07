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

SAMPLES = ["TTV", "VV", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]

#DIR_MC= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC_FullRazorInclusive/combined"
DIR_MC= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC"
DIR_DATA= "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data"
PREFIX = "RazorInclusive"
DATA_NAME='SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_razorskim'
FILENAMES = {
        "TTJets"    : DIR_MC+"/"+PREFIX+"_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root",
        "WJets"     : DIR_MC+"/"+PREFIX+"_WJetsToLNu_HTBinned_1pb_weighted.root",
        #"WJets"     : DIR_MC+"/"+PREFIX+"_WJets_1pb_weighted.root",
        "SingleTop" : DIR_MC+"/"+PREFIX+"_SingleTop_1pb_weighted.root",
        "VV" : DIR_MC+"/"+PREFIX+"_VV_1pb_weighted.root",
        "TTV" : DIR_MC+"/"+PREFIX+"_TTV_1pb_weighted.root",
        "DYJets"     : DIR_MC+"/"+PREFIX+"_DYJetsToLL_M-50_HTBinned_1pb_weighted.root",
        #"DYJets"     : DIR_MC+"/"+PREFIX+"_DYJetsToLL_1pb_weighted.root",
        "ZInv"     : DIR_MC+"/"+PREFIX+"_ZJetsToNuNu_HTBinned_1pb_weighted.root",
        #"ZInv"     : DIR_MC+"/"+PREFIX+"_ZJetsToNuNu_1pb_weighted.root",
        "Data"      : DIR_DATA+'/'+PREFIX+'_'+DATA_NAME+'.root',
        }

config = "config/run2_sideband.config"
FIT_DIR = "MyFitResults"
#config = "config/backup.config"
#FIT_DIR = "FitResultsSideband"
TOYS_FILES = {
        "MuMultiJet":FIT_DIR+"/toys_Bayes_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/toys_Bayes_EleMultiJet.root",
        }

weightOpts = ["doNVtxWeights"]
#shapeErrors = ["muoneff", "eleeff", "jes"]
#miscErrors = ["mt"]
shapeErrors = []
miscErrors = []

cfg = Config.Config(config)
binsMRLep = cfg.getBinning("MuMultiJet")[0]
binsRsqLep = cfg.getBinning("MuMultiJet")[1]
leptonicBinning = { "MR":binsMRLep, "Rsq":binsRsqLep }
blindBins = [(x,y) for x in range(3,len(leptonicBinning["MR"])+1) for y in range(2,len(leptonicBinning["Rsq"])+1)]

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

    #get scale factor histograms
    sfHists = loadScaleFactorHists(processNames=SAMPLES, debugLevel=debugLevel)

    #estimate yields in signal region
    for lepType in ["Mu", "Ele"]:
        for jets in ["MultiJet"]:
            boxName = lepType+jets
            btaglist = [0]
            #btaglist = [0,1,2,3]
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
                        bins=leptonicSignalRegionBins, lumiMC=MCLUMI, lumiData=LUMI, 
                        weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                        weightOpts=weightOpts, shapeErrors=shapeErrors, miscErrors=miscErrors,
                        fitToyFiles=TOYS_FILES, boxName=boxName, blindBins=blindBins,
                        btags=nBtags, debugLevel=debugLevel)
