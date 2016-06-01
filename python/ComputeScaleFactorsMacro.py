import sys
import argparse
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import wjetsSingleLeptonCutsMC, ttjetsSingleLeptonCutsMC, xbinsSignal, colsSignal, cutsData_2016, cutsData_2015
from macro.razorWeights import *
from macro.razorMacros import *
from SidebandMacro import LUMI as LUMI_DATA

MCLUMI = 1 

SAMPLES_TTJ1L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["Other", "QCD", "DYJets", "SingleTop", "TTJets", "WJets"]
SAMPLES_TTJ2L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_DYJ2L = ["Other", "SingleTop", "WJets", "TTJets", "DYJets"]
SAMPLES_Veto = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L_INV = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "TTJets", "WJetsInv"]
SAMPLES_DYJ2L_INV = ["Other", "SingleTop", "WJets", "TTJets", "DYJets"]

DIR_1L = "ControlRegions"
#DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonFull_1p23_2015Final/RazorSkim/"
PREFIX_1L = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
FILENAMES_1L = {
            "TTJets"   : DIR_1L+"/"+PREFIX_1L+"_TTJets_1pb_weighted_RazorSkim.root",
            "WJets"    : DIR_1L+"/"+PREFIX_1L+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_1L+"/"+PREFIX_1L+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_1L+"/"+PREFIX_1L+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_1L+"/"+PREFIX_1L+"_Other_1pb_weighted_RazorSkim.root",
            "ZInv"     : DIR_1L+"/"+PREFIX_1L+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
            "QCD"      : DIR_1L+"/"+PREFIX_1L+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",       
            "Data"     : DIR_1L+"/"+PREFIX_1L+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }

DIR_2L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFull_1p23_2015Final/RazorSkim/"
PREFIX_2L = "RunTwoRazorControlRegions_DileptonFull_DileptonSkim"
FILENAMES_2L = {
            "TTJets"   : DIR_2L+"/"+PREFIX_2L+"_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root",
            "WJets"    : DIR_2L+"/"+PREFIX_2L+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_2L+"/"+PREFIX_2L+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_2L+"/"+PREFIX_2L+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_2L+"/"+PREFIX_2L+"_VV_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_2L+"/"+PREFIX_2L+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }

DIR_1L_INV = "ControlRegions"
#DIR_1L_INV = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonAddToMETFull_1p23_2015Final/RazorSkim/"
PREFIX_1L_INV = "RunTwoRazorControlRegions_OneLeptonAddToMetFull_SingleLeptonSkim"
FILENAMES_1L_INV = {
            "TTJets"   : DIR_1L_INV+"/"+PREFIX_1L_INV+"_TTJets_1pb_weighted_RazorSkim.root",
            "WJetsInv" : DIR_1L_INV+"/"+PREFIX_1L_INV+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_1L_INV+"/"+PREFIX_1L_INV+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_1L_INV+"/"+PREFIX_1L_INV+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_1L_INV+"/"+PREFIX_1L_INV+"_Other_1pb_weighted_RazorSkim.root",
            "ZInv"     : DIR_1L_INV+"/"+PREFIX_1L_INV+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
            "QCD"      : DIR_1L_INV+"/"+PREFIX_1L_INV+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_1L_INV+"/"+PREFIX_1L_INV+"_SingleLepton_Run2015D_RazorSkim_GoodLumiGolden_NoDuplicates.root"
            }

DATA_2016 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p1/OneLeptonFull/RunTwoRazorControlRegions_OneLeptonFull_SingleLepton_2016B_PRv2_GoodLumiGolden_RazorSkim_NoDuplicates.root"

weightOpts = []

config = "config/run2_20151229_ControlRegion.config"
cfg = Config.Config(config)
binsMRLep = cfg.getBinning("WJetControlRegion")[0]
binsRsqLep = cfg.getBinning("WJetControlRegion")[1]
binsNBTags = [0.,1.,2.,3.,4.]
binsNJets80 = [0.,1.,2.,20.]
binsNJets = [2.,3.,4.,5.,6.,7.,20.]
binsLepPt = [20.,25.,30.,35.,40.,45.,50.,70.,100]
ControlRegionBinning = { "MR":binsMRLep, "Rsq":binsRsqLep, "NBJetsMedium":binsNBTags, "NJets80":binsNJets80, "NJets40":binsNJets, "lep1.Pt()":binsLepPt, ("MR","Rsq"):[], ("MR","Rsq","NBJetsMedium"):[]}
ZNuNu_1L_ControlRegionBinning = { "MR_NoW":binsMRLep, "Rsq_NoW":binsRsqLep, "NBJetsMedium":binsNBTags, "NJets80":binsNJets80, "NJets40":binsNJets, ("MR_NoW","Rsq_NoW"):[] }

xbinsWJETS1L = xbinsSignal["WJetControlRegion"]["0B"]
colsWJETS1L = colsSignal["WJetControlRegion"]["0B"]
xbinsTTJETS1L = xbinsSignal["TTJetsSingleLeptonControlRegion"]["0B"]
colsTTJETS1L = colsSignal["TTJetsSingleLeptonControlRegion"]["0B"]
xbinsWJETS1LINV = xbinsSignal["WJetInvControlRegion"]["0B"]
colsWJETS1LINV = colsSignal["WJetInvControlRegion"]["0B"]

printdir="ControlSamplePlots"

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--use-2016-data", dest="use2016Data", action="store_true",
                                help="Use 2016 data")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug

    #initialize
    weightHists = {}
    sfHists = {}
    cutsData = cutsData_2015

    #choose data file
    if args.use2016Data: 
        FILENAMES_1L["Data"] = DATA_2016
        FILENAMES_2L["Data"] = DATA_2016
        FILENAMES_1L_INV["Data"] = DATA_2016
        printdir = "ControlSamplePlots2016"
        LUMI_DATA = 589
        cutsData = cutsData_2016

    #make output directory
    os.system('mkdir -p '+printdir)

    sfVars = ("MR","Rsq")
    sfVars_NoW = ("MR_NoW", "Rsq_NoW")

    plotOpts = { 'comment':False }

    ##########################################################
    #TTJets control sample
    ##########################################################
    ttjetsSingleLeptonHists = makeControlSampleHists("TTJetsSingleLepton", 
                filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
                cutsMC=ttjetsSingleLeptonCutsMC, cutsData=cutsData["TTJetsSingleLepton"], 
                bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
                weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, 
                printdir=printdir, sfVars=sfVars, debugLevel=debugLevel,
                unrollBins=(xbinsTTJETS1L, colsTTJETS1L), plotOpts=plotOpts)
    appendScaleFactors("TTJets", ttjetsSingleLeptonHists, sfHists, lumiData=LUMI_DATA, th2PolyXBins=xbinsTTJETS1L, th2PolyCols=colsTTJETS1L, debugLevel=debugLevel, var=sfVars, printdir=printdir)
    macro.exportHists(ttjetsSingleLeptonHists, outFileName='controlHistogramsTTJetsSingleLepton.root', outDir=printdir, debugLevel=debugLevel)

    ##########################################################
    #WJets control sample
    ##########################################################
    #sfHists_ForWJets = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorScaleFactors_Inclusive_TTJets.root", processNames=SAMPLES_WJ1L, debugLevel=debugLevel)
    wjetsSingleLeptonHists = makeControlSampleHists("WJetsSingleLepton", 
                filenames=FILENAMES_1L, samples=SAMPLES_WJ1L, 
                cutsMC=wjetsSingleLeptonCutsMC, cutsData=cutsData["WJetsSingleLepton"], 
                bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
                weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, 
                printdir=printdir, plotDensity=False, sfVars=sfVars, debugLevel=debugLevel,
                unrollBins=(xbinsWJETS1L, colsWJETS1L), plotOpts=plotOpts)
    appendScaleFactors("WJets", wjetsSingleLeptonHists, sfHists, lumiData=LUMI_DATA, th2PolyXBins=xbinsWJETS1L, th2PolyCols=colsWJETS1L, debugLevel=debugLevel, var=sfVars, printdir=printdir)
    macro.exportHists(wjetsSingleLeptonHists, outFileName='controlHistogramsWJetsSingleLepton.root', outDir=printdir, debugLevel=debugLevel)


    ##########################################################
    #W+Jets Add lepton to MET control sample
    ##########################################################    
    #sfHists_ForWJetsInv = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorScaleFactors_Inclusive_TTJets.root", processNames=SAMPLES_WJ1L_INV, debugLevel=debugLevel)
    wjetsSingleLeptonInvHists = makeControlSampleHists("WJetsSingleLeptonInv", 
                 filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
                 cutsMC=wjetsSingleLeptonInvCutsMC, cutsData=cutsData["WJetsSingleLeptonInv"], 
                 bins=ZNuNu_1L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
                 weightHists=weightHists, sfHists=sfHists, weightOpts=weightOpts, 
                 printdir=printdir, plotDensity=False, sfVars=sfVars, debugLevel=debugLevel,
                 unrollBins=(xbinsWJETS1LINV, colsWJETS1LINV), plotOpts=plotOpts)
    appendScaleFactors("WJetsInv", wjetsSingleLeptonInvHists, sfHists, var=sfVars_NoW, lumiData=LUMI_DATA, th2PolyXBins=xbinsWJETS1LINV, th2PolyCols=colsWJETS1LINV, debugLevel=debugLevel, printdir=printdir)
    macro.exportHists(wjetsSingleLeptonInvHists, outFileName='controlHistogramsWJetsSingleLeptonInv.root', outDir=printdir, debugLevel=debugLevel)

    #write scale factors
    outfile = rt.TFile("RazorScaleFactors.root", "RECREATE")
    for name in sfHists:
        print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
        sfHists[name].Write(sfHists[name].GetName().replace("Poly",""))
    outfile.Close()
