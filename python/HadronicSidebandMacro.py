import sys,os
import argparse
import ROOT as rt

#local imports
from framework import Config
import macro.macro as macro
from macro.razorAnalysis import *
from macro.razorWeights import *
from macro.razorMacros import *

LUMI = 1264 #in /pb
MCLUMI = 1 

DIR_DATA="root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data"
PREFIX = "RazorInclusive"
DATA_NAME='HTMHT_Run2015D_Oct05ReMiniAOD_PRv4_GoodLumiGolden'
FILENAMES = {
        "Data"      : DIR_DATA+'/'+PREFIX+'_'+DATA_NAME+'.root',
        }

config = "config/run2_sideband.config"
FIT_DIR = "MyFitResults"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/toys_Bayes_MultiJet.root",
        }

cfg = Config.Config(config) 
binsMR = cfg.getBinning("MultiJet")[0]
binsRsq = cfg.getBinning("MultiJet")[1]
hadronicBinning = { "MR":binsMR, "Rsq":binsRsq }

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

    #estimate yields in signal region
    boxName = "MultiJet"
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
                filenames=FILENAMES, samples=[], 
                cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                bins=hadronicBinning, lumiMC=MCLUMI, lumiData=LUMI, 
                weightHists={}, sfHists={}, treeName="RazorInclusive", 
                weightOpts=[], shapeErrors=[], miscErrors=[],
                fitToyFiles=TOYS_FILES, boxName=boxName, blindBins=hadronicBlindBins,
                btags=nBtags, debugLevel=debugLevel)
