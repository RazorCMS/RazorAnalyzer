import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
import macro.macro as macro
from macro.razorAnalysis import razorCuts
from macro.razorWeights import loadScaleFactorHists 
from macro.razorMacros import runFitAndToys, makeControlSampleHists

LUMI = 2185 #in /pb
MCLUMI = 1 

SAMPLES = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"]
BOXES = ["MuMultiJet", "MultiJet", "EleMultiJet"]

DIR_MC = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108/"
DIR_DATA = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106"
DATA_NAMES={
    'MultiJet':'RazorInclusive_HTMHT_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered',
    'EleMultiJet':'RazorInclusive_SingleElectron_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered',
    'MuMultiJet':'RazorInclusive_SingleMuon_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered',
    }
FILENAMES_MC = {
        "TTJets"    : DIR_MC+"/"+"FullRazorInclusive_TTJets_1pb_weighted.root",
        "WJets"     : DIR_MC+"/"+"FullRazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop" : DIR_MC+"/"+"FullRazorInclusive_SingleTop_1pb_weighted.root",
        "Other" : DIR_MC+"/"+"FullRazorInclusive_Other_1pb_weighted.root",
        "DYJets"     : DIR_MC+"/"+"FullRazorInclusive_DYJetsToLL_M-50_HTBinned_1pb_weighted.root",
        "ZInv"     : DIR_MC+"/"+"FullRazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root",
        }
FILENAMES = {name:copy.copy(FILENAMES_MC) for name in BOXES}
for name in BOXES: FILENAMES[name]["Data"] = DIR_DATA+'/'+DATA_NAMES[name]+'.root'

config = "config/run2_20151108_Preapproval_2b3b_data.config"
FIT_DIR = "eos/cms/store/group/phys_susy/razor/Run2Analysis/FitResults/ResultForDecemberJamboree2015/Data_2093ipb"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/MultiJet/sideband/toys_Bayes_MultiJet.root",
        "MuMultiJet":FIT_DIR+"/MuMultiJet/sideband/toys_Bayes_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/EleMultiJet/sideband/toys_Bayes_EleMultiJet.root",
        }

weightOpts = []
shapeErrors = []
miscErrors = []

cfg = Config.Config(config)
binsMRHad = cfg.getBinning("MultiJet")[0]
binsRsqHad = cfg.getBinning("MultiJet")[1]
hadronicBinning = { "MR":binsMRHad, "Rsq":binsRsqHad, ("MR","Rsq"):[] }
binsMRLep = cfg.getBinning("MuMultiJet")[0]
binsRsqLep = cfg.getBinning("MuMultiJet")[1]
leptonicBinning = { "MR":binsMRLep, "Rsq":binsRsqLep, ("MR","Rsq"):[] }
binning = { "MultiJet":hadronicBinning, "MuMultiJet":leptonicBinning, "EleMultiJet":leptonicBinning}

blindBins = {b:[(x,y) for x in range(2,len(binning[b]["MR"])+1) for y in range(2,len(binning[b]["Rsq"])+1)] for b in binning}

dirName="SignalRegionPlots"

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--unblind", help="do not blind signal sensitive region", action='store_true')
    parser.add_argument('--no-mc', help="do not process MC, do data and fit only", action='store_true', dest="noMC")
    parser.add_argument('--no-fit', help="do not load fit results, process data and MC only", action='store_true', dest='noFit')
    parser.add_argument('--full', help="do full fit (default is sideband)", action='store_true')
    parser.add_argument('--no-data', help="do not process data, do fit and MC only", action='store_true', dest='noData')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug

    plotOpts = {}
    doSideband=(not args.full)
    if not doSideband:
        FIT_DIR = FIT_DIR.replace('Sideband','Full').replace('sideband','full')
        TOYS_FILES = {b:TOYS_FILES[b].replace('sideband','full').replace('Sideband','full') for b in TOYS_FILES}
        dirName += '_Full'
        plotOpts['sideband'] = False
    else:
        plotOpts['sideband'] = True
    if args.unblind:
        dirName += '_Unblinded'
    if args.noFit: 
        TOYS_FILES = None
        del plotOpts['sideband']

    #initialize
    weightHists = {}
    sfHists = {}

    #make output directory
    os.system('mkdir -p '+dirName)

    #get scale factor histograms
    sfHists = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root", processNames=SAMPLES, scaleFactorNames={"ZInv":"WJetsInv"}, debugLevel=debugLevel)

    #estimate yields in signal region
    for boxName in BOXES:

        #apply options
        blindBinsToUse = blindBins[boxName]
        if args.unblind: blindBinsToUse = None
        samplesToUse = SAMPLES
        if args.noMC: samplesToUse = []
        if samplesToUse is None or len(samplesToUse) == 0:
            filesToUse = {"Data":FILENAMES[boxName]["Data"]}
        else:
            filesToUse = FILENAMES[boxName]
        if args.noData: 
            del filesToUse['Data']

        #loop over btag bins
        btaglist = [0]
        #btaglist = [0,1,2,3]
        for btags in btaglist:
            #if btags == 3: continue #temporary
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
            if not args.noFit and not os.path.isfile(TOYS_FILES[boxName]):
                print "Fit file",TOYS_FILES[boxName],"not found, trying to recreate it"
                runFitAndToys(FIT_DIR, boxName, LUMI, DATA_NAMES[boxName], DIR_DATA, config=config, sideband=doSideband)
                #check
                if not os.path.isfile(TOYS_FILES[boxName]):
                    print "Error creating fit file",TOYS_FILES[boxName]
                    sys.exit()
            makeControlSampleHists(extboxName, 
                    filenames=filesToUse, samples=samplesToUse, 
                    cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                    bins=binning[boxName], lumiMC=MCLUMI, lumiData=LUMI, 
                    weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                    weightOpts=weightOpts, shapeErrors=shapeErrors, miscErrors=miscErrors,
                    fitToyFiles=TOYS_FILES, boxName=boxName, blindBins=blindBinsToUse,
                    btags=nBtags, debugLevel=debugLevel, printdir=dirName, plotOpts=plotOpts)
