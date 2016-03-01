import sys
import argparse
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import wjetsSingleLeptonCutsMC, wjetsSingleLeptonCutsData, ttjetsSingleLeptonCutsMC, ttjetsSingleLeptonCutsData
from macro.razorWeights import *
from macro.razorMacros import *
from ComputeScaleFactorsMacro import xbinsWJETS1L, colsWJETS1L, xbinsTTJETS1L, colsTTJETS1L, xbinsWJETS1LINV, colsWJETS1LINV
from SidebandMacro import LUMI as LUMI_DATA

MCLUMI = 1 

SAMPLES_TTJ1L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "TTJets", "WJets"]
SAMPLES_TTJ2L = ["Other", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_DYJ2L = ["Other", "SingleTop", "WJets", "TTJets", "DYJets"]
SAMPLES_Veto = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"]
SAMPLES_WJ1L_INV = ["Other", "ZInv", "QCD", "DYJets", "SingleTop", "TTJets", "WJetsInv"]
SAMPLES_DYJ2L_INV = ["Other", "SingleTop", "WJets", "TTJets", "DYJets"]
#apply correct version of MR when applying scale factors
ScaleFactorVars_WJETS1L_INV = { 
                              "WJetsInv":("MR_NoW","Rsq_NoW"),
                              "TTJets":("MR","Rsq"),
                              }

DIR_1L = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonFull_1p23_2015Final/RazorSkim/"
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

DIR_1L_INV = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/OneLeptonAddToMETFull_1p23_2015Final/RazorNJets80Skim/"
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

DIR_2L_INV = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/DileptonFullAddToMET_1p23_2015Final/RazorNJets80Skim/"
PREFIX_2L_INV = "RunTwoRazorControlRegions_DileptonAddToMetFull_DileptonSkim"
FILENAMES_2L_INV = {
            "TTJets"   : DIR_2L_INV+"/"+PREFIX_2L_INV+"_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root",
            #"WJets"    : DIR_2L_INV+"/"+PREFIX_2L_INV+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
            "SingleTop": DIR_2L_INV+"/"+PREFIX_2L_INV+"_SingleTop_1pb_weighted_RazorSkim.root",
            "DYJets"   : DIR_2L_INV+"/"+PREFIX_2L_INV+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
            "Other"    : DIR_2L_INV+"/"+PREFIX_2L_INV+"_Other_1pb_weighted_RazorSkim.root",
            "Data"     : DIR_2L_INV+"/"+PREFIX_2L_INV+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
            }


weightOpts = []

config = "config/run2_20151229_ControlRegion.config"
cfg = Config.Config(config)
binsMRLep = cfg.getBinning("WJetControlRegion")[0]
binsRsqLep = cfg.getBinning("WJetControlRegion")[1]
binsNBTags = [0.,1.,2.,3.,4.]
binsNJets80 = [0.,1.,2.,9.]
binsNJets = [0.,1.,2.,3.,4.,7.,20.]
binsLepPt = [20.,25.,30.,35.,40.,45.,50.,70.,100]
ControlRegionBinning = { "MR":binsMRLep, "Rsq":binsRsqLep, "NBJetsMedium":binsNBTags, "NJets80":binsNJets80, "NJets40":binsNJets, "lep1.Pt()":binsLepPt, ("MR","Rsq"):[]}
ZNuNu_1L_ControlRegionBinning = { "MR_NoW":binsMRLep, "Rsq_NoW":binsRsqLep, "NBJetsMedium":binsNBTags, "NJets80":binsNJets80, "NJets40":binsNJets, ("MR_NoW","Rsq_NoW"):[] }

binsMR_7Jet = [400,700,900,4000]
binsRsq_7Jet = [0.15,0.25,0.30,1.5]
ControlRegionBinning7Jet = { "MR":binsMR_7Jet, "Rsq":binsRsq_7Jet, "NBJetsMedium":binsNBTags, ("MR","Rsq"):[] }

printdir="TestScaleFactorPlots"

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
    sfHists = {}

    #make output directory
    os.system('mkdir -p '+printdir)

    sfVars = ("MR","Rsq")
    sfVars_NoW = ("MR_NoW", "Rsq_NoW")

    plotOpts = { "comment":False }

    ##########################################################
    # #DYJets control sample
    ##########################################################
    # dyjetsDileptonHists = makeControlSampleHists("DYJetsDilepton", 
    #             filenames=FILENAMES_2L, samples=SAMPLES_DYJ2L, 
    #             cutsMC=dyjetsDileptonCutsMC, cutsData=dyjetsDileptonCutsData, 
    #             bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #             weightHists=weightHists, plotDensity=False, sfHists=sfHists, weightOpts=weightOpts, 
    #             printdir=printdir, debugLevel=debugLevel)
    # appendScaleFactors("DYJets", dyjetsDileptonHists, sfHists, lumiData=LUMI_DATA, debugLevel=debugLevel, printdir=printdir)

    ##########################################################
    #1L W+Jets & TTJets control sample
    ##########################################################    
    #sfHists_OneLeptonScaleFactorClosureTest_Uncorr = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root", processNames=SAMPLES_TTJ1L, debugLevel=debugLevel)
    #OneLeptonScaleFactorClosureTestHists_Uncorr = makeControlSampleHists("OneLeptonScaleFactorClosureTestUncorr", 
    #              filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
    #              cutsMC=OneLeptonScaleFactorClosureTestMC, cutsData=OneLeptonScaleFactorClosureTestData, 
    #              bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #              weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest_Uncorr, weightOpts=weightOpts, 
    #              printdir=printdir, btags=-1, plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
    #              plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
    #sfHists_OneLeptonScaleFactorClosureTest = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root", processNames=SAMPLES_TTJ1L, debugLevel=debugLevel)
    #OneLeptonScaleFactorClosureTestHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest", 
    #              filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
    #              cutsMC=OneLeptonScaleFactorClosureTestMC, cutsData=OneLeptonScaleFactorClosureTestData, 
    #              bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #              weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest, weightOpts=weightOpts, 
    #              printdir=printdir, btags=-1, plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
    #              plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
    #OneLeptonScaleFactorClosureTest0BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_0B", 
    #           filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
    #           cutsMC=OneLeptonScaleFactorClosureTest0BMC, cutsData=OneLeptonScaleFactorClosureTest0BData, 
    #           bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest, weightOpts=weightOpts, 
    #           printdir=printdir, btags=-1, plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
    #           plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
    #appendScaleFactors("OneLepton0BMR", OneLeptonScaleFactorClosureTest0BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton0BRsq", OneLeptonScaleFactorClosureTest0BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)
    #OneLeptonScaleFactorClosureTest1BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_1B", 
    #           filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
    #           cutsMC=OneLeptonScaleFactorClosureTest1BMC, cutsData=OneLeptonScaleFactorClosureTest1BData, 
    #           bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest, weightOpts=weightOpts, 
    #           printdir=printdir, btags=-1, plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
    #           plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
    #appendScaleFactors("OneLepton1BMR", OneLeptonScaleFactorClosureTest1BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton1BRsq", OneLeptonScaleFactorClosureTest1BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)
    #OneLeptonScaleFactorClosureTest2BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_2B", 
    #           filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
    #           cutsMC=OneLeptonScaleFactorClosureTest2BMC, cutsData=OneLeptonScaleFactorClosureTest2BData, 
    #           bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest, weightOpts=weightOpts, 
    #           printdir=printdir, btags=-1, plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
    #           plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
    #appendScaleFactors("OneLepton2BMR", OneLeptonScaleFactorClosureTest2BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton2BRsq", OneLeptonScaleFactorClosureTest2BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)
    #OneLeptonScaleFactorClosureTest3BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_3B", 
    #           filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
    #           cutsMC=OneLeptonScaleFactorClosureTest3BMC, cutsData=OneLeptonScaleFactorClosureTest3BData, 
    #           bins=ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest, weightOpts=weightOpts, 
    #           printdir=printdir, btags=-1, plotDensity=True, sfVars=sfVars, debugLevel=debugLevel,
    #           plotOpts=plotOpts, unrollBins=(xbinsTTJETS1L, colsTTJETS1L))
    #appendScaleFactors("OneLepton3BMR", OneLeptonScaleFactorClosureTest3BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton3BRsq", OneLeptonScaleFactorClosureTest3BHists, sfHists_OneLeptonScaleFactorClosureTest, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)
  
     #Cross Check for >= 7 Jets bin
   # sfHists_OneLeptonScaleFactorClosureTest7Jet = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedTo7Jet.root", processNames=SAMPLES_TTJ1L, debugLevel=debugLevel)
   # OneLeptonScaleFactorClosureTest7JetHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_7JetBin", 
   #               filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
   #               cutsMC=OneLeptonScaleFactorClosureTest7JetMC, cutsData=OneLeptonScaleFactorClosureTest7JetData, 
   #               bins=ControlRegionBinning7Jet, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
   #               weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest7Jet, weightOpts=weightOpts, 
   #               printdir=printdir, btags=-1, plotDensity=False, sfVars=sfVars, debugLevel=debugLevel)

   # OneLeptonScaleFactorClosureTest7Jet0BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_7JetBin0B", 
   #               filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
   #               cutsMC=OneLeptonScaleFactorClosureTest7Jet0BMC, cutsData=OneLeptonScaleFactorClosureTest7Jet0BData, 
   #               bins=ControlRegionBinning7Jet, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
   #               weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest7Jet, weightOpts=weightOpts, 
   #               printdir=printdir, btags=-1, plotDensity=False, sfVars=sfVars, debugLevel=debugLevel)
    #appendScaleFactors("OneLepton7Jet0BMR", OneLeptonScaleFactorClosureTest7Jet0BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton7Jet0BRsq", OneLeptonScaleFactorClosureTest7Jet0BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)

   # OneLeptonScaleFactorClosureTest7Jet1BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_7JetBin1B", 
   #               filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
   #               cutsMC=OneLeptonScaleFactorClosureTest7Jet1BMC, cutsData=OneLeptonScaleFactorClosureTest7Jet1BData, 
   #               bins=ControlRegionBinning7Jet, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
   #               weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest7Jet, weightOpts=weightOpts, 
   #               printdir=printdir, btags=-1, plotDensity=False, sfVars=sfVars, debugLevel=debugLevel)
    #appendScaleFactors("OneLepton7Jet1BMR", OneLeptonScaleFactorClosureTest7Jet1BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton7Jet1BRsq", OneLeptonScaleFactorClosureTest7Jet1BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)

   # OneLeptonScaleFactorClosureTest7Jet2BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_7JetBin2B", 
   #               filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
   #               cutsMC=OneLeptonScaleFactorClosureTest7Jet2BMC, cutsData=OneLeptonScaleFactorClosureTest7Jet2BData, 
   #               bins=ControlRegionBinning7Jet, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
   #               weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest7Jet, weightOpts=weightOpts, 
   #               printdir=printdir, btags=-1, plotDensity=False, sfVars=sfVars, debugLevel=debugLevel)
    #appendScaleFactors("OneLepton7Jet2BMR", OneLeptonScaleFactorClosureTest7Jet2BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton7Jet2BRsq", OneLeptonScaleFactorClosureTest7Jet2BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)

   # OneLeptonScaleFactorClosureTest7Jet3BHists = makeControlSampleHists("OneLeptonScaleFactorClosureTest_7JetBin3B", 
   #               filenames=FILENAMES_1L, samples=SAMPLES_TTJ1L, 
   #               cutsMC=OneLeptonScaleFactorClosureTest7Jet3BMC, cutsData=OneLeptonScaleFactorClosureTest7Jet3BData, 
   #               bins=ControlRegionBinning7Jet, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
   #               weightHists=weightHists, sfHists=sfHists_OneLeptonScaleFactorClosureTest7Jet, weightOpts=weightOpts, 
   #               printdir=printdir, btags=-1, plotDensity=False, sfVars=sfVars, debugLevel=debugLevel)
    #appendScaleFactors("OneLepton7Jet3BMR", OneLeptonScaleFactorClosureTest7Jet3BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLepton7Jet3BRsq", OneLeptonScaleFactorClosureTest7Jet3BHists, sfHists_OneLeptonScaleFactorClosureTest7Jet, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq", signifThreshold=1.0, printdir=printdir)

    #outfile = rt.TFile("Razor7JetBTagClosureTests.root", "RECREATE")
    #for name in sfHists_OneLeptonScaleFactorClosureTest7Jet:
    #    if "OneLepton7Jet" in name:
    #        print "Writing scale factor histogram",sfHists_OneLeptonScaleFactorClosureTest7Jet[name].GetName(),"to file"
    #        sfHists_OneLeptonScaleFactorClosureTest7Jet[name].Write()
    #outfile.Close()


    # ##########################################################
    # #W+Jets Add lepton to MET control sample
    # ##########################################################
    sfHists_wjetsInv = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_Uncorrected.root", processNames=SAMPLES_WJ1L_INV, debugLevel=debugLevel)
    wjetsSingleLeptonInvHists = makeControlSampleHists("OneLeptonInvScaleFactorClosureTestUncorr", 
                filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
                cutsMC=OneLeptonInvScaleFactorClosureTestMC, cutsData=OneLeptonInvScaleFactorClosureTestData, 
                bins=ZNuNu_1L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
                weightHists=weightHists, sfHists=sfHists_wjetsInv, weightOpts=weightOpts, 
                printdir=printdir, plotDensity=True, sfVars=ScaleFactorVars_WJETS1L_INV , debugLevel=debugLevel)
    #sfHists_wjetsInv = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root", processNames=SAMPLES_WJ1L_INV, debugLevel=debugLevel)
    #wjetsSingleLeptonInvHists = makeControlSampleHists("OneLeptonInvScaleFactorClosureTest", 
    #            filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
    #            cutsMC=OneLeptonInvScaleFactorClosureTestMC, cutsData=OneLeptonInvScaleFactorClosureTestData, 
    #            bins=ZNuNu_1L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #            weightHists=weightHists, sfHists=sfHists_wjetsInv, weightOpts=weightOpts, 
    #            printdir=printdir, plotDensity=True, sfVars=ScaleFactorVars_WJETS1L_INV , debugLevel=debugLevel)
    #wjetsSingleLeptonInvHists = makeControlSampleHists("OneLeptonInvScaleFactorClosureTest0B", 
    #             filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
    #             cutsMC=OneLeptonInvScaleFactorClosureTest0BMC, cutsData=OneLeptonInvScaleFactorClosureTest0BData, 
    #             bins=ZNuNu_1L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #             weightHists=weightHists, sfHists=sfHists_wjetsInv, weightOpts=weightOpts, 
    #             printdir=printdir, plotDensity=True, sfVars=ScaleFactorVars_WJETS1L_INV , debugLevel=debugLevel,
    #             plotOpts=plotOpts, unrollBins=(xbinsWJETS1LINV, colsWJETS1LINV))
    #appendScaleFactors("OneLeptonInv0BMR", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR_NoW", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLeptonInv0BRsq", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq_NoW", signifThreshold=1.0, printdir=printdir)
    #wjetsSingleLeptonInvHists = makeControlSampleHists("OneLeptonInvScaleFactorClosureTest1B", 
    #           filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
    #           cutsMC=OneLeptonInvScaleFactorClosureTest1BMC, cutsData=OneLeptonInvScaleFactorClosureTest1BData, 
    #           bins=ZNuNu_1L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, sfHists=sfHists_wjetsInv, weightOpts=weightOpts, 
    #           printdir=printdir, plotDensity=True, sfVars=ScaleFactorVars_WJETS1L_INV , debugLevel=debugLevel,
    #           plotOpts=plotOpts, unrollBins=(xbinsWJETS1LINV, colsWJETS1LINV))
    #appendScaleFactors("OneLeptonInv1BMR", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR_NoW", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLeptonInv1BRsq", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq_NoW", signifThreshold=1.0, printdir=printdir)
    #wjetsSingleLeptonInvHists = makeControlSampleHists("OneLeptonInvScaleFactorClosureTest2B", 
    #           filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
    #           cutsMC=OneLeptonInvScaleFactorClosureTest2BMC, cutsData=OneLeptonInvScaleFactorClosureTest2BData, 
    #           bins=ZNuNu_1L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, sfHists=sfHists_wjetsInv, weightOpts=weightOpts, 
    #           printdir=printdir, plotDensity=True, sfVars=ScaleFactorVars_WJETS1L_INV , debugLevel=debugLevel,
    #           plotOpts=plotOpts, unrollBins=(xbinsWJETS1LINV, colsWJETS1LINV))
    #appendScaleFactors("OneLeptonInv2BMR", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR_NoW", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLeptonInv2BRsq", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq_NoW", signifThreshold=1.0, printdir=printdir)
    #wjetsSingleLeptonInvHists = makeControlSampleHists("OneLeptonInvScaleFactorClosureTest3B", 
    #           filenames=FILENAMES_1L_INV, samples=SAMPLES_WJ1L_INV, 
    #           cutsMC=OneLeptonInvScaleFactorClosureTest3BMC, cutsData=OneLeptonInvScaleFactorClosureTest3BData, 
    #           bins=ZNuNu_1L_ControlRegionBinning, lumiMC=MCLUMI, lumiData=LUMI_DATA, 
    #           weightHists=weightHists, sfHists=sfHists_wjetsInv, weightOpts=weightOpts, 
    #           printdir=printdir, plotDensity=True, sfVars=ScaleFactorVars_WJETS1L_INV , debugLevel=debugLevel,
    #           plotOpts=plotOpts, unrollBins=(xbinsWJETS1LINV, colsWJETS1LINV))
    #appendScaleFactors("OneLeptonInv3BMR", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="MR_NoW", signifThreshold=1.0, printdir=printdir)
    #appendScaleFactors("OneLeptonInv3BRsq", wjetsSingleLeptonInvHists, sfHists_wjetsInv, lumiData=LUMI_DATA, debugLevel=debugLevel, var="Rsq_NoW", signifThreshold=1.0, printdir=printdir)


    ##write scale factors
    #outfile = rt.TFile("RazorBTagClosureTests.root", "RECREATE")
    #for name in sfHists_OneLeptonScaleFactorClosureTest:
    #    if "OneLepton" in name:
    #        print "Writing scale factor histogram",sfHists_OneLeptonScaleFactorClosureTest[name].GetName(),"to file"
    #        sfHists_OneLeptonScaleFactorClosureTest[name].Write()
    #for name in sfHists_wjetsInv:
    #    if "OneLepton" in name:
    #        print "Writing scale factor histogram",sfHists_wjetsInv[name].GetName(),"to file"
    #        sfHists_wjetsInv[name].Write()
    #outfile.Close()

 
