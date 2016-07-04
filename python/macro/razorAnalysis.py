## Inclusive razor analysis definitions
import ROOT as rt
from array import array
import sys
import copy
import csv

#local imports
from framework import Config

#####################################
### WEIGHTS NEEDED
#####################################

razorWeightOpts = {
        "Razor2015":[],
        "Razor2016":[
                     #"reapplyNPUWeights", #remove pileup weight and multiply by new PU weight
                     #"nbjets" #apply an extra btag correction
                     ], 
        }
razorWeightHists = {
        "Razor2015":{},
        "Razor2016":{ "pileup":
            ("data/PileupWeights/PileupReweight2016_06172016.root", "PileupReweight")
                }
        }

#####################################
### NTUPLES
#####################################

razorNtuples = { 
        "SingleLepton":{},
        "SingleLeptonInv":{},
        "DileptonInv":{},
        "VetoLepton":{},
        "VetoTau":{},
        "SignalHadronic":{},
        "SignalMuon":{},
        "SignalElectron":{},
        "SignalLepton":{},
        }

razorSamples = {
        "TTJetsSingleLepton":["Other", "DYJets", "SingleTop", "WJets", "TTJets"],
        "WJetsSingleLepton":["Other", "DYJets", "SingleTop", "TTJets", "WJets"],
        #"WJetsSingleLepton":["Other", "QCD", "DYJets", "SingleTop", "TTJets", "WJets"],
        "TTJetsDilepton":["Other", "DYJets", "SingleTop", "WJets", "TTJets"],
        "VetoLepton":["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"],
        "VetoTau":["Other", "ZInv", "QCD", "DYJets", "SingleTop", "WJets", "TTJets"],
        "WJetsSingleLeptonInv":["DYJets", "SingleTop", "TTJets", "WJetsInv"],
        #"WJetsSingleLeptonInv":["Other", "ZInv", "QCD", "DYJets", "SingleTop", "TTJets", "WJetsInv"],
        "DYJetsDileptonInv":["Other", "SingleTop", "WJets", "TTJets", "DYJetsInv"],
        "SignalHadronic":["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"],
        "SignalLeptonic":["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets"],
        #"SignalHadronic":["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets2L", "TTJets1L"],
        #"SignalLeptonic":["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"],
        }
razorSamplesReduced = {
        "TTJetsSingleLepton":["Other", "WJets", "TTJets"],
        "WJetsSingleLepton":["Other", "TTJets", "WJets"],
        #"WJetsSingleLepton":["Other", "QCD", "TTJets", "WJets"],
        "TTJetsDilepton":["Other", "WJets", "TTJets"],
        "VetoLepton":["Other", "ZInv", "QCD", "WJets", "TTJets"],
        "VetoTau":["Other", "ZInv", "QCD", "WJets", "TTJets"],
        "WJetsSingleLeptonInv":["TTJets", "WJetsInv"],
        #"WJetsSingleLeptonInv":["Other", "ZInv", "QCD", "TTJets", "WJetsInv"],
        "DYJetsDileptonInv":["Other", "WJets", "TTJets", "DYJetsInv"],
        "SignalHadronic":["Other", "QCD", "ZInv", "WJets", "TTJets"],
        "SignalLeptonic":["Other", "ZInv", "WJets", "TTJets"],
        }

### 2015 ntuples
dir1L2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/OneLeptonFull_1p23_2015Final/RazorSkim/"
dir1LInv2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/OneLeptonAddToMETFull_1p23_2015Final/RazorSkim/"
dir2LInv2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/DileptonFullAddToMET_1p23_2015Final/RazorNJets80Skim/"
dirVetoL2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/VetoLeptonFull_1p23_2015Final/RazorNJets80Skim/"
dirVetoTau2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2015/VetoTauFull_1p23_2015Final/"
dirSignalMC2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2015/V1p23_Background_20160108/"
dirSignalData2015 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/2015/V1p23_ForMoriond20160119/RazorSkim/"

prefix1L2015 = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
prefix1LInv2015 = "RunTwoRazorControlRegions_OneLeptonAddToMetFull_SingleLeptonSkim" 
prefix2LInv2015 = "RunTwoRazorControlRegions_DileptonAddToMetFull_DileptonSkim"
prefixVetoL2015 = "RunTwoRazorControlRegions_VetoLeptonFull"
prefixVetoTau2015 = "RunTwoRazorControlRegions_VetoTauFull_RazorSkim"

### 1-lepton control region
razorNtuples["SingleLepton"]["Razor2015"] = {
        "TTJets"   : dir1L2015+"/"+prefix1L2015+"_TTJets_1pb_weighted_RazorSkim.root",
        "WJets"    : dir1L2015+"/"+prefix1L2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop": dir1L2015+"/"+prefix1L2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "DYJets"   : dir1L2015+"/"+prefix1L2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "Other"    : dir1L2015+"/"+prefix1L2015+"_Other_1pb_weighted_RazorSkim.root",
        "ZInv"     : dir1L2015+"/"+prefix1L2015+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
        "QCD"      : dir1L2015+"/"+prefix1L2015+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",       
        "Data"     : dir1L2015+"/"+prefix1L2015+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
        }
 
### 1-lepton invisible control region
razorNtuples["SingleLeptonInv"]["Razor2015"] = {
        "TTJets"   : dir1LInv2015+"/"+prefix1LInv2015+"_TTJets_1pb_weighted_RazorSkim.root",
        "WJetsInv" : dir1LInv2015+"/"+prefix1LInv2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop": dir1LInv2015+"/"+prefix1LInv2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "DYJets"   : dir1LInv2015+"/"+prefix1LInv2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "Other"    : dir1LInv2015+"/"+prefix1LInv2015+"_Other_1pb_weighted_RazorSkim.root",
        "ZInv"     : dir1LInv2015+"/"+prefix1LInv2015+"_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
        "QCD"      : dir1LInv2015+"/"+prefix1LInv2015+"_QCD_HTBinned_1pb_weighted_RazorSkim.root",
        "Data"     : dir1LInv2015+"/"+prefix1LInv2015+"_SingleLepton_Run2015D_RazorSkim_GoodLumiGolden_NoDuplicates.root"
        }

### 2-lepton invisible control region
razorNtuples["DileptonInv"]["Razor2015"] = {
        "TTJets"    : dir2LInv2015+"/"+prefix2LInv2015+"_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted_RazorSkim.root",
        "WJets"     : dir2LInv2015+"/"+prefix2LInv2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop" : dir2LInv2015+"/"+prefix2LInv2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "DYJetsInv" : dir2LInv2015+"/"+prefix2LInv2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
        "Other"     : dir2LInv2015+"/"+prefix2LInv2015+"_Other_1pb_weighted_RazorSkim.root",
        "Data"      : dir2LInv2015+"/"+prefix2LInv2015+"_SingleLepton_Run2015D_GoodLumiGolden_NoDuplicates_RazorSkim.root"
        }

### Veto lepton control region
razorNtuples["VetoLepton"]["Razor2015"] = {
        "TTJets"   : dirVetoL2015+"/"+prefixVetoL2015+"_TTJets_1pb_weighted_RazorSkim.root",
        "WJets"    : dirVetoL2015+"/"+prefixVetoL2015+"_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
        "SingleTop": dirVetoL2015+"/"+prefixVetoL2015+"_SingleTop_1pb_weighted_RazorSkim.root",
        "Data"     : dirVetoL2015+"/"+prefixVetoL2015+"_HTMHT_Run2015D_GoodLumiGolden_RazorSkim.root"
        }

### Veto tau control region
razorNtuples["VetoTau"]["Razor2015"] = {
        "TTJets"   : dirVetoTau2015+"/"+prefixVetoTau2015+"_TTJets_1pb_weighted.root",
        "WJets"    : dirVetoTau2015+"/"+prefixVetoTau2015+"_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop": dirVetoTau2015+"/"+prefixVetoTau2015+"_SingleTop_1pb_weighted.root",
        "DYJets"   : dirVetoTau2015+"/"+prefixVetoTau2015+"_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted.root",
        "Other"    : dirVetoTau2015+"/"+prefixVetoTau2015+"_Other_1pb_weighted.root",
        "Data"     : dirVetoTau2015+"/"+prefixVetoTau2015+"_HTMHT_Run2015D_GoodLumiGolden.root"
        }

### Signal region
razorNtuples["SignalHadronic"]["Razor2015"] = {
        "TTJets1L" : dirSignalMC2015+"/FullRazorInclusive_TTJets1L_1pb_weighted.root",
        "TTJets2L" : dirSignalMC2015+"/FullRazorInclusive_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
        "WJets"    : dirSignalMC2015+"/FullRazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop": dirSignalMC2015+"/FullRazorInclusive_SingleTop_1pb_weighted.root",
        "Other"    : dirSignalMC2015+"/FullRazorInclusive_Other_1pb_weighted.root",
        "DYJets"   : dirSignalMC2015+"/FullRazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted.root",
        "ZInv"     : dirSignalMC2015+"/FullRazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root",
        #data-driven QCD prediction for MultiJet
        "QCD"      : dirSignalData2015+"/RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root", 
        "Data"     : dirSignalData2015+"/RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"
        }

razorNtuples["SignalMuon"]["Razor2015"] = razorNtuples["SignalHadronic"]["Razor2015"].copy()

razorNtuples["SignalMuon"]["Razor2015"]["Data"] = dirSignalData2015+"RazorInclusive_SingleMuon_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"

razorNtuples["SignalElectron"]["Razor2015"] = razorNtuples["SignalHadronic"]["Razor2015"].copy()

razorNtuples["SignalElectron"]["Razor2015"]["Data"] = dirSignalData2015+"RazorInclusive_SingleElectron_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"

razorNtuples["SignalLepton"]["Razor2015"] = razorNtuples["SignalHadronic"]["Razor2015"].copy()

razorNtuples["SignalLepton"]["Razor2015"]["Data"] = dirSignalData2015+"RazorInclusive_SingleLepton_Run2015D_GoodLumiGolden_RazorSkim_CSCBadTrackFilter.root"

### 2016 ntuples
dirCR2016 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/"
dirSR2016 = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/"
versionMC2016 = "V3p2"
versionData2016 = "V3p3"

prefix1L2016 = "RunTwoRazorControlRegions_OneLeptonFull_SingleLeptonSkim"
prefix1LInv2016 = "RunTwoRazorControlRegions_OneLeptonAddToMetFull_SingleLeptonSkim" 
prefix2LInv2016 = "RunTwoRazorControlRegions_DileptonAddToMetFull_DileptonSkim"
prefixVetoL2016 = "RunTwoRazorControlRegions_VetoLeptonFull"
prefixVetoTau2016 = "RunTwoRazorControlRegions_VetoTauFull_RazorSkim"
skimstr = "_RazorSkim"

#on EOS
#dir1L2016 = dirCR2016+'/'+versionMC2016+'/OneLeptonRazorSkimLeptonSkim'
#dir1LInv2016 = dirCR2016+'/'+versionMC2016+'/OneLeptonInvRazorSkimLeptonSkim'
#dir2LInv2016 = dirCR2016+'/'+versionMC2016+'/DileptonRazorSkimDileptonSkim'
#dirVetoL2016 = dirCR2016+'/'+versionMC2016+'/VetoLeptonRazorSkim'
#dirVetoTau2016 = dirCR2016+'/'+versionMC2016+'/VetoTauRazorSkim'
#dirSignal2016 = dirSR2016+'/'+versionMC2016

#local directories
dir1L2016 = 'Backgrounds/1L'
dir1LInv2016 = 'Backgrounds/1LInv'
dir2LInv2016 = 'Backgrounds/2LInv'
dirVetoL2016 = 'Backgrounds/VetoL'
dirVetoTau2016 = 'Backgrounds/VetoTau'
dirSignal2016 = 'Backgrounds/Signal'

razorNtuples["SingleLepton"]["Razor2016"] = {
        "TTJets"   : dir1L2016+"/"+prefix1L2016+"_TTJets_1pb_weighted"+skimstr+".root",
        "WJets"    : dir1L2016+"/"+prefix1L2016+"_WJets_1pb_weighted"+skimstr+".root",
        "SingleTop": dir1L2016+"/"+prefix1L2016+"_SingleTop_1pb_weighted"+skimstr+".root",
        "DYJets"   : dir1L2016+"/"+prefix1L2016+"_DYJets_1pb_weighted"+skimstr+".root",
        "Other"    : dir1L2016+"/"+prefix1L2016+"_Other_1pb_weighted"+skimstr+".root",
        "Data"     : dir1L2016+"/"+prefix1L2016+"_Data_NoDuplicates_GoodLumiGolden.root"
        }
razorNtuples["SingleLeptonInv"]["Razor2016"] = {
        "TTJets"   : dir1LInv2016+"/"+prefix1LInv2016+"_TTJets_1pb_weighted"+skimstr+".root",
        "WJetsInv"    : dir1LInv2016+"/"+prefix1LInv2016+"_WJets_1pb_weighted"+skimstr+".root",
        "SingleTop": dir1LInv2016+"/"+prefix1LInv2016+"_SingleTop_1pb_weighted"+skimstr+".root",
        "DYJets"   : dir1LInv2016+"/"+prefix1LInv2016+"_DYJets_1pb_weighted"+skimstr+".root",
        #"Other"    : dir1LInv2016+"/"+prefix1LInv2016+"_Other_1pb_weighted"+skimstr+".root",
        "Data"     : dir1LInv2016+"/"+prefix1LInv2016+"_Data_NoDuplicates_GoodLumiGolden.root"
        }
razorNtuples["DileptonInv"]["Razor2016"] = {
        "TTJets"   : dir2LInv2016+"/"+prefix2LInv2016+"_TTJets_1pb_weighted"+skimstr+".root",
        "WJets"    : dir2LInv2016+"/"+prefix2LInv2016+"_WJets_1pb_weighted"+skimstr+".root",
        "SingleTop": dir2LInv2016+"/"+prefix2LInv2016+"_SingleTop_1pb_weighted"+skimstr+".root",
        "DYJetsInv"   : dir2LInv2016+"/"+prefix2LInv2016+"_DYJets_1pb_weighted"+skimstr+".root",
        "Other"    : dir2LInv2016+"/"+prefix2LInv2016+"_Other_1pb_weighted"+skimstr+".root",
        #"ZInv"     : dir2LInv2016+"/"+prefix2LInv2016+"_ZInv_1pb_weighted"+skimstr+".root",
        #"QCD"      : dir2LInv2016+"/"+prefix2LInv2016+"_QCD_1pb_weighted"+skimstr+".root",       
        "Data"     : dir2LInv2016+"/"+prefix2LInv2016+"_Data_NoDuplicates_GoodLumiGolden.root"
        }
razorNtuples["VetoLepton"]["Razor2016"] = {
        "TTJets"   : dirVetoL2016+"/"+prefixVetoL2016+"_TTJets_1pb_weighted"+skimstr+".root",
        "WJets"    : dirVetoL2016+"/"+prefixVetoL2016+"_WJets_1pb_weighted"+skimstr+".root",
        "SingleTop": dirVetoL2016+"/"+prefixVetoL2016+"_SingleTop_1pb_weighted"+skimstr+".root",
        "DYJets"   : dirVetoL2016+"/"+prefixVetoL2016+"_DYJets_1pb_weighted"+skimstr+".root",
        "Other"    : dirVetoL2016+"/"+prefixVetoL2016+"_Other_1pb_weighted"+skimstr+".root",
        "ZInv"     : dirVetoL2016+"/"+prefixVetoL2016+"_ZInv_1pb_weighted"+skimstr+".root",
        "QCD"      : dirVetoL2016+"/"+prefixVetoL2016+"_QCD_1pb_weighted"+skimstr+".root",       
        "Data"     : dirVetoL2016+"/"+prefixVetoL2016+"_Data_NoDuplicates_GoodLumiGolden.root"
        }
razorNtuples["VetoTau"]["Razor2016"] = {
        "TTJets"   : dirVetoTau2016+"/"+prefixVetoTau2016+"_TTJets_1pb_weighted"+skimstr+".root",
        "WJets"    : dirVetoTau2016+"/"+prefixVetoTau2016+"_WJets_1pb_weighted"+skimstr+".root",
        "SingleTop": dirVetoTau2016+"/"+prefixVetoTau2016+"_SingleTop_1pb_weighted"+skimstr+".root",
        "DYJets"   : dirVetoTau2016+"/"+prefixVetoTau2016+"_DYJets_1pb_weighted"+skimstr+".root",
        "Other"    : dirVetoTau2016+"/"+prefixVetoTau2016+"_Other_1pb_weighted"+skimstr+".root",
        "ZInv"     : dirVetoTau2016+"/"+prefixVetoTau2016+"_ZInv_1pb_weighted"+skimstr+".root",
        "QCD"      : dirVetoTau2016+"/"+prefixVetoTau2016+"_QCD_1pb_weighted"+skimstr+".root",       
        "Data"     : dirVetoTau2016+"/"+prefixVetoTau2016+"_Data_NoDuplicates_GoodLumiGolden.root"
        }
razorNtuples["SignalHadronic"]["Razor2016"] = {
        "TTJets"   : dirSignal2016+"/FullRazorInclusive_TTJets_1pb_weighted"+skimstr+".root",
        "WJets"    : dirSignal2016+"/FullRazorInclusive_WJets_1pb_weighted"+skimstr+".root",
        "SingleTop": dirSignal2016+"/FullRazorInclusive_SingleTop_1pb_weighted"+skimstr+".root",
        "DYJets"   : dirSignal2016+"/FullRazorInclusive_DYJets_1pb_weighted"+skimstr+".root",
        "Other"    : dirSignal2016+"/FullRazorInclusive_Other_1pb_weighted"+skimstr+".root",
        "ZInv"     : dirSignal2016+"/FullRazorInclusive_ZInv_1pb_weighted"+skimstr+".root",
        #QCD predicted using data driven method
        "QCD"      : dirSignal2016+"/FullRazorInclusive_Data_NoDuplicates_GoodLumiGolden2p6fb.root",
        "Data"     : dirSignal2016+"/FullRazorInclusive_Data_NoDuplicates_GoodLumiGolden2p6fb.root"
        #"QCD"      : dirSignal2016+"/FullRazorInclusive_Data_NoDuplicates_GoodLumiGolden800pb.root",
        #"Data"     : dirSignal2016+"/FullRazorInclusive_Data_NoDuplicates_GoodLumiGolden800pb.root"
        }
razorNtuples["SignalLepton"]["Razor2016"] = razorNtuples["SignalHadronic"]["Razor2016"].copy()
razorNtuples["SignalMuon"]["Razor2016"] = razorNtuples["SignalHadronic"]["Razor2016"].copy()
razorNtuples["SignalElectron"]["Razor2016"] = razorNtuples["SignalHadronic"]["Razor2016"].copy()

#update version number for data, in case different
for name,files in razorNtuples.iteritems():
    if "Razor2016" in files:
        files["Razor2016"]["Data"] = files["Razor2016"]["Data"].replace(
                '/'+versionMC2016+'/','/'+versionData2016+'/')

#####################################
### TRIGGER
#####################################

triggerNames = {}
for ttype in ["Razor", "SingleLepton", "Dilepton"]: triggerNames[ttype] = { "Data":{}, "MC":{} }

# 2015 triggers
# Data
triggerNames["Razor"]["Data"]["Razor2015"] = [
        'HLT_RsqMR240_Rsq0p09_MR200',
        'HLT_RsqMR240_Rsq0p09_MR200_4jet',
        'HLT_RsqMR270_Rsq0p09_MR200',
        'HLT_RsqMR270_Rsq0p09_MR200_4jet',
        'HLT_RsqMR260_Rsq0p09_MR200',
        'HLT_RsqMR260_Rsq0p09_MR200_4jet',
        'HLT_RsqMR300_Rsq0p09_MR200',
        'HLT_RsqMR300_Rsq0p09_MR200_4jet',
        'HLT_Rsq0p25',
        'HLT_Rsq0p30',
        'HLT_Rsq0p36',
        ]
triggerNames["SingleLepton"]["Data"]["Razor2015"] = [
        'HLT_Mu50',
        'HLT_IsoMu20',
        'HLT_IsoMu27',
        'HLT_IsoTkMu20',
        'HLT_IsoTkMu27',
        'HLT_Ele23_WPLoose_Gsf',
        'HLT_Ele27_WPLoose_Gsf',
        'HLT_Ele27_eta2p1_WPLoose_Gsf',
        'HLT_Ele27_eta2p1_WPTight_Gsf',
        'HLT_Ele32_eta2p1_WPLoose_Gsf',
        'HLT_Ele32_eta2p1_WPTight_Gsf',
        'HLT_Ele105_CaloIdVT_GsfTrkIdT',
        'HLT_Ele115_CaloIdVT_GsfTrkIdT',
        ]
triggerNames["Dilepton"]["Data"]["Razor2015"] = [
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ',
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
        'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
        'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
        'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL',
        ]
# MC
triggerNames["Razor"]["MC"]["Razor2015"] = copy.copy(triggerNames["Razor"]["Data"]["Razor2015"])
triggerNames["SingleLepton"]["MC"]["Razor2015"] = [
        'HLT_Mu50',
        'HLT_IsoMu20',
        'HLT_IsoMu27',
        'HLT_IsoTkMu20',
        'HLT_IsoTkMu27',
        'HLT_Ele23_WP75_Gsf',
        'HLT_Ele27_WP85_Gsf',
        'HLT_Ele27_eta2p1_WP75_Gsf',
        'HLT_Ele32_eta2p1_WP75_Gsf',
        'HLT_Ele105_CaloIdVT_GsfTrkIdT',
        'HLT_Ele115_CaloIdVT_GsfTrkIdT',
        'HLT_Ele22_eta2p1_WP75_Gsf',
        ]
triggerNames["Dilepton"]["MC"]["Razor2015"] = copy.copy(triggerNames["Dilepton"]["Data"]["Razor2015"])

# 2016 trigger
# Data
triggerNames["Razor"]["Data"]["Razor2016"] = copy.copy(triggerNames["Razor"]["Data"]["Razor2015"])
triggerNames["SingleLepton"]["Data"]["Razor2016"] = [
        'HLT_Mu50',
        'HLT_Mu45_eta2p1',
        'HLT_IsoMu22',
        'HLT_IsoTkMu22',
        'HLT_Ele25_eta2p1_WPTight_Gsf',
        'HLT_Ele27_WPTight_Gsf',
        'HLT_Ele27_eta2p1_WPLoose_Gsf',
        'HLT_Ele45_WPLoose_Gsf',
        'HLT_Ele105_CaloIdVT_GsfTrkIdT',
        ]
triggerNames["Dilepton"]["Data"]["Razor2016"] = copy.copy(triggerNames["Dilepton"]["Data"]["Razor2015"])
#MC (samples have broken HLT menu)
triggerNames["Razor"]["MC"]["Razor2016"] = [ 'PASS' ]
triggerNames["SingleLepton"]["MC"]["Razor2016"] = [ 'PASS' ]
triggerNames["Dilepton"]["MC"]["Razor2016"] = [ 'PASS' ]

class TriggerUtils(object):
    """Stores and retrieves trigger names"""
    def __init__(self, trigType, tag, isData=False):
        self.trigType = trigType
        self.tag = tag
        if isData:
            self.isDataTag = "Data"
        else:
            self.isDataTag = "MC"

        # get trigger path name/number correspondence
        self.triggerDict = {}
        if tag == 'Razor2015':
            trigFile = 'data/RazorHLTPathnames_2015.dat'
        elif tag == 'Razor2016':
            trigFile = 'data/RazorHLTPathnames_2016.dat'
        else:
            self.badInit(trigType,tag)
        with open(trigFile) as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                self.triggerDict[row[-1]] = int(row[0]) 

        # get trigger names
        if self.trigType in triggerNames and self.isDataTag in triggerNames[self.trigType] and self.tag in triggerNames[self.trigType][self.isDataTag]:
            self.names = triggerNames[self.trigType][self.isDataTag][self.tag]
        else:
            self.badInit(trigType,tag)

        # get trigger numbers
        self.getTriggerNums()

    def badInit(self, trigType, tag):
        print "Error in triggerUtils: trigger type/tag combination '",trigType,tag,"' not recognized!"
        self.names = None

    def getTriggerNums(self):
        self.nums = []
        for name in self.names:
            if name == 'PASS': # indicates that no trigger cuts should be applied
                self.nums = [-1] 
                return
            elif name in self.triggerDict:
                self.nums.append( self.triggerDict[name] )
            else:
                print "Warning in triggerUtils.getTriggerNums: trigger name",name,"not found!"

    def appendTriggerCuts(self, cutsString):
        """Append a string of the form "(HLTDecision[t1] || HLTDecision[t2] || ... || HLTDecision[tN]) && " 
           to the provided cut string, where t1...tN are the desired trigger numbers"""
        if -1 in self.nums: # trigger requirement not applied
            return cutsString
        return '('+(' || '.join(['HLTDecision['+str(n)+']' for n in self.nums]))+") && "+cutsString

#####################################
### CUTS
#####################################

razorCuts = {}

def appendBoxCuts(cuts, boxNums):
    """Append a string of the form "(box == b1 || box == b2 || ... || box == bN) && " to the provided cut string, where b1...bN are the desired box numbers"""
    return '('+(' || '.join(['box == '+str(n) for n in boxNums])) + ") && " + cuts

recommendedNoiseFilters = ["Flag_HBHENoiseFilter","Flag_HBHEIsoNoiseFilter","Flag_goodVertices",
        "Flag_eeBadScFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_CSCTightHaloFilter",
        #"Flag_badChargedCandidateFilter","Flag_badMuonFilter"
        ]
def appendNoiseFilters(cuts, tree=None):
    ret = copy.copy(cuts)
    for bit in recommendedNoiseFilters:
        if (tree is None) or hasattr(tree, bit):
            ret += " && " + bit
    return ret

# 1-lepton control region
razorCuts["SingleLepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 25) || (abs(lep1Type) == 13 && lep1.Pt() > 20)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && MR > 300 && Rsq > 0.15"

# 1-lepton invisible control region
razorCuts["SingleLeptonInv"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && ((abs(lep1Type) == 11 && lep1.Pt() > 25) || (abs(lep1Type) == 13 && lep1.Pt() > 20)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && NJets80_NoW >= 2 && MR_NoW > 300 && Rsq_NoW > 0.15 && ( weight < 0.01 || weight == 1)"

# TTJets 2-lepton control region
razorCuts["TTJetsDilepton"] = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1PassTight && lep2PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 20)) && ((abs(lep2Type) == 11 && lep2.Pt() > 25) || (abs(lep2Type) == 13 && lep2.Pt() > 20)) && mll > 20 && ((abs(lep1Type) != abs(lep2Type)) || (mll < 76 || mll > 106)) && NBJetsMedium > 0 && MET > 40 && MR > 300 && Rsq > 0.15"

# DYJets 2-lepton invisible control region
razorCuts["DYJetsDileptonInv"] = "((abs(lep1Type) == 11 && abs(lep2Type) == 11) || (abs(lep1Type) == 13 && abs(lep2Type) == 13)) && lep1.Pt() > 30 && lep2.Pt() > 20 && recoZmass > 80 && recoZmass < 110 && NBJetsMedium == 0 && NJets80_NoZ >= 2 && MR_NoZ > 300 && Rsq_NoZ > 0.15"

# Veto lepton control region
razorCuts["VetoLepton"]  = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

razorCuts["VetoElectron"] = "(abs(lep1Type) == 11) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

razorCuts["VetoMuon"] = "(abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

razorCuts["VetoTau"] = "(abs(lep1Type) == 15) && lep1PassLoose && lep1.Pt() > 20 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 400 && Rsq > 0.25"

# Hadronic boxes with inverted lepton veto
cutsMultiJetForVeto = "(box == 11 || box == 12 || box == 21) && MR > 400.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"
cutsDiJetForVeto = "(box == 14) && MR > 400.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

razorCuts["MultiJetForVetoLepton"] = cutsMultiJetForVeto+" && ( abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13 ) && leadingGenLeptonPt > 5"

razorCuts["MultiJetForVetoTau"] = cutsMultiJetForVeto+" && abs(leadingGenLeptonType) == 15 && leadingGenLeptonPt > 20"

razorCuts["DiJetForVetoLepton"] = cutsDiJetForVeto+" && ( abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13 ) && leadingGenLeptonPt > 5"

razorCuts["DiJetForVetoTau"] = cutsDiJetForVeto+" && abs(leadingGenLeptonType) == 15 && leadingGenLeptonPt > 20"

# Signal boxes
razorBoxes = {
        "MuEle" : [0], 
        "MuMu" : [1],
        "EleEle" : [2],
        "MuSixJet" : [3],
        "MuFourJet" : [4],
        "MuJet" : [5],
        "EleSixJet" : [6],
        "EleFourJet" : [7],
        "EleJet" : [8],
        "LeptonSixJet" : [3,6],
        "LeptonFourJet" : [4,7],
        "LeptonJet" : [5,8],
        "LooseLeptonSixJet" : [9],
        "LooseLeptonFourJet" : [10],
        "LooseLeptonDiJet" : [13],
        "SixJet" : [11],
        "FourJet" : [12],
        "DiJet" : [14],	  
        "MuMultiJet" : [3,4,18],
        "EleMultiJet" : [6,7,19],
        "LeptonMultiJet" : [3,4,18,6,7,19],
        "LooseLeptonMultiJet" : [9,10,20],
        "MultiJet" : [11,12,21],
        "FourToSixJet" : [11,12,21],
        "SevenJet" : [11,12,21],
        }

hadronicRazorBoxes = ["DiJet", "FourJet", "SixJet", "MultiJet", "FourToSixJet", "SevenJet"]
looseLeptonRazorBoxes = ["LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet", "LooseLeptonMultiJet"]
leptonicRazorBoxes = ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet","EleJet", "EleFourJet", "EleSixJet", "EleMultiJet", "LeptonJet", "LeptonFourJet", "LeptonSixJet", "LeptonMultiJet"]
electronRazorBoxes = ["EleJet", "EleFourJet", "EleSixJet", "EleMultiJet"]
muonRazorBoxes = ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet"]
leptonRazorBoxes = ["LeptonJet", "LeptonFourJet", "LeptonSixJet", "LeptonMultiJet"]
dileptonRazorBoxes = ["MuEle", "MuMu", "EleEle"]

dileptonSignalRegionCuts = "MR > 400.000000 && MR < 4000 && Rsq > 0.150000 && Rsq < 1.5 && abs(dPhiRazor) < 2.8"
leptonicSignalRegionCuts = "MR > 400.000000 && MR < 4000 && Rsq > 0.150000 && Rsq < 1.5 && mT > 120"
looseLeptonSignalRegionCuts = "MR > 500.000000 && MR < 4000 && Rsq > 0.250000 && Rsq < 1.5 && mTLoose > 100 && nJets80 >= 2"
hadronicSignalRegionCuts = "MR > 500.000000 && MR < 4000 && Rsq > 0.250000 && Rsq < 1.5 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

for box in razorBoxes:
    if box in hadronicRazorBoxes: 
        razorCuts[box] = appendBoxCuts(hadronicSignalRegionCuts, razorBoxes[box])
    elif box in looseLeptonRazorBoxes:
        razorCuts[box] = appendBoxCuts(looseLeptonSignalRegionCuts, razorBoxes[box])
    elif box in leptonicRazorBoxes:
        razorCuts[box] = appendBoxCuts(leptonicSignalRegionCuts, razorBoxes[box])
    elif box in dileptonRazorBoxes:
        razorCuts[box] = appendBoxCuts(dileptonSignalRegionCuts, razorBoxes[box])
razorCuts['FourToSixJet'] += ' && nSelectedJets < 7'
razorCuts['SevenJet'] += ' && nSelectedJets >= 7'

#####################################
### BINNING
#####################################

razorBinning = {}

### Single Lepton Control Region
razorBinning["SingleLepton"] = {
        "MR" : [300, 400, 500, 600, 700, 900, 1200, 4000],
        "Rsq" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets40" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        }
razorBinning["SingleLeptonInv"] = {
        "MR_NoW" : [300, 400, 500, 600, 700, 900, 1200, 4000],
        "Rsq_NoW" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets_NoW" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        }

### Veto lepton control region
razorBinning["VetoLepton"] =  {
        "MR" : [400, 500, 600, 700, 900, 1200, 4000],
        "Rsq" : [0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets40" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "lep1.Pt()" : [5, 10, 15, 20.,30.,40.,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }
razorBinning["SignalRegionForVetoLepton"] =  {
        "MR" : [400, 500, 600, 700, 900, 1200, 4000],
        "Rsq" : [0.25,0.30,0.41,0.52,0.64,1.5],
        "nSelectedJets" : [1,2,4,20],
        "leadingGenLeptonPt" : [5, 10, 15, 20.,30.,40.,100,1000],
        "abs(leadingGenLeptonEta)" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }

### Veto Tau Control Region
razorBinning["VetoTau"] = {
        "MR" : [400, 500, 600, 700, 900, 4000],
        "Rsq" : [0.25,0.30,0.41,1.5],
        "NJets40" : [1,2,4,20],
        "NBJetsMedium" : [0,1,2,3,4],
        "lep1.Pt()" : [20,30,40,100,1000],
        "abs(lep1.Eta())" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }
razorBinning["SignalRegionForVetoTau"] = {
        "MR" : [400, 500, 600, 700, 900, 4000],
        "Rsq" : [0.25,0.30,0.41,1.5],
        "nSelectedJets" : [1,2,4,20],
        "leadingGenLeptonPt" : [20,30,40,100,1000],
        "abs(leadingGenLeptonEta)" : [0, 0.5, 1.0, 1.5, 2.0, 2.5],
        }

### DYJets Dilepton Invisible Control Region
razorBinning["DYJetsDileptonInv"] = {
        "MR_NoZ" : [400, 500, 600, 700, 900, 1200, 4000],
        "Rsq_NoZ" : [0.25,0.30,0.41,0.52,0.64,1.5],
        "NJets_NoZ" : [1,2,4,20]
        }
razorBinning["DYJetsDileptonInvReduced"] = {
        "MR_NoZ" : [400, 500, 600, 700, 900, 4000],
        "Rsq_NoZ" : [0.25,0.30,0.41,1.5],
        "NJets_NoZ" : [1,2,4,20]
        }

### Signal region
signalConfig = "config/run2_20151108_Preapproval_2b3b_data.config"
cfg = Config.Config(signalConfig)
razorBinning["SignalHadronic"] = {
        "MR" : cfg.getBinning("MultiJet")[0],
        "Rsq" : cfg.getBinning("MultiJet")[1],
        #"leadingJetPt": range(0, 1000, 100),
        #"subleadingJetPt": range(0, 1000, 100),
        }
razorBinning["SignalLeptonic"] = {
        "MR" : cfg.getBinning("MuMultiJet")[0],
        "Rsq" : cfg.getBinning("MuMultiJet")[1],
        #"leadingJetPt": range(0, 1000, 100),
        #"subleadingJetPt": range(0, 1000, 100),
        }

### Add 2D razor binning to each region
for region,binning in razorBinning.iteritems():
    if "MR" in binning and "Rsq" in binning:
        binning[("MR","Rsq")] = []
    if "MR_NoW" in binning and "Rsq_NoW" in binning:
        binning[("MR_NoW","Rsq_NoW")] = []
    if "MR_NoZ" in binning and "Rsq_NoZ" in binning:
        binning[("MR_NoZ","Rsq_NoZ")] = []
    if "MR_NoPho" in binning and "Rsq_NoPho" in binning:
        binning[("MR_NoPho","Rsq_NoPho")] = []

### Non-grid binning 
xbinsSignal = { "MultiJet":{}, "MuMultiJet":{}, "EleMultiJet":{}, "LeptonMultiJet":{},
                "DiJet":{}, "MuJet":{}, "EleJet":{}, "FourToSixJet":{}, "SevenJet":{}}
colsSignal = { "MultiJet":{}, "MuMultiJet":{}, "EleMultiJet":{}, "LeptonMultiJet":{},
                "DiJet":{}, "MuJet":{}, "EleJet":{}, "FourToSixJet":{}, "SevenJet":{}}

xbinsSignal["MultiJet"]["0B"] = [ 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MultiJet"]["1B"] = [ 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MultiJet"]["2B"] = [ 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MultiJet"]["3B"] = [ 500, 600, 700, 900, 1200, 1600, 4000 ]

xbinsSignal["FourToSixJet"] = xbinsSignal["MultiJet"]

xbinsSignal["MuMultiJet"]["0B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MuMultiJet"]["1B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MuMultiJet"]["2B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["MuMultiJet"]["3B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]


xbinsSignal["EleMultiJet"]["0B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["EleMultiJet"]["1B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["EleMultiJet"]["2B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]
xbinsSignal["EleMultiJet"]["3B"] = [ 400, 500, 600, 700, 900, 1200, 1600, 4000 ]


xbinsSignal["DiJet"] = xbinsSignal["MultiJet"]
xbinsSignal["MuJet"] = xbinsSignal["MuMultiJet"]
xbinsSignal["EleJet"] = xbinsSignal["EleMultiJet"]
xbinsSignal["LeptonMultiJet"] = xbinsSignal["MuMultiJet"]
xbinsSignal["LeptonJet"] = xbinsSignal["LeptonMultiJet"]

colsSignal["MultiJet"]["0B"] = [
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        ]
colsSignal["MultiJet"]["1B"] = [
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        ]
colsSignal["MultiJet"]["2B"] = [
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 1.5 ],
        ]
colsSignal["MultiJet"]["3B"] = [
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 0.52, 0.64, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        ]
colsSignal["MuMultiJet"]["0B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        ]
colsSignal["MuMultiJet"]["1B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        ]
colsSignal["MuMultiJet"]["2B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        ]
colsSignal["MuMultiJet"]["3B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        ]
colsSignal["EleMultiJet"]["0B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        ]
colsSignal["EleMultiJet"]["1B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        ]
colsSignal["EleMultiJet"]["2B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        ]
colsSignal["EleMultiJet"]["3B"] = [
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 0.30, 0.41, 0.52, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15,0.20, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        [ 0.15, 0.25, 1.5 ],
        ]

colsSignal["FourToSixJet"] = colsSignal["MultiJet"]

colsSignal["DiJet"] = colsSignal["MultiJet"]
colsSignal["MuJet"] = colsSignal["MuMultiJet"]
colsSignal["EleJet"] = colsSignal["EleMultiJet"]
colsSignal["LeptonMultiJet"] = colsSignal["MuMultiJet"]
colsSignal["LeptonJet"] = colsSignal["LeptonMultiJet"]

xbinsSignal["SevenJet"]["0B"] = [ 500, 600, 700, 900, 1200, 4000 ]
xbinsSignal["SevenJet"]["1B"] = [ 500, 600, 700, 900, 1200, 4000 ]
xbinsSignal["SevenJet"]["2B"] = [ 500, 600, 700, 900, 1200, 4000 ]
xbinsSignal["SevenJet"]["3B"] = [ 500, 600, 700, 4000 ]

colsSignal["SevenJet"]["0B"] = [
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        ]
colsSignal["SevenJet"]["1B"] = [
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        ]
colsSignal["SevenJet"]["2B"] = [
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 1.5 ],
        ]
colsSignal["SevenJet"]["3B"] = [
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        [ 0.25, 0.30, 0.41, 1.5 ],
        ]

#binning for scale factor histograms
xbinsSignal["TTJetsSingleLepton"] = {}
colsSignal["TTJetsSingleLepton"] = {}
xbinsSignal["TTJetsSingleLepton"]["0B"] = [300, 400, 500, 600, 700, 900, 1200, 4000]
colsSignal["TTJetsSingleLepton"]["0B"] = [
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 1.5 ],
    [ 0.15, 0.2, 0.25, 1.5 ],
    [ 0.15, 0.2, 1.5 ],
    ]

xbinsSignal["WJetsSingleLepton"] = {}
colsSignal["WJetsSingleLepton"] = {}
xbinsSignal["WJetsSingleLepton"]["0B"] = [300, 400, 500, 600, 700, 900, 1200, 4000]
colsSignal["WJetsSingleLepton"]["0B"] = [
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 1.5 ],
    [ 0.15, 0.2, 0.25, 1.5 ],
    ]

xbinsSignal["WJetsSingleLeptonInv"] = {}
colsSignal["WJetsSingleLeptonInv"] = {}
xbinsSignal["WJetsSingleLeptonInv"]["0B"] = [300, 400, 500, 600, 700, 900, 1200, 4000]
colsSignal["WJetsSingleLeptonInv"]["0B"] = [
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.3, 0.41, 1.5 ],
    [ 0.15, 0.2, 0.25, 0.30, 0.41, 1.5 ],
    ]

xbinsSignal["GJetsInv"] = {}
colsSignal["GJetsInv"] = {}
xbinsSignal["GJetsInv"]["0B"] = [400, 500, 600, 700, 900, 1200, 4000]
colsSignal["GJetsInv"]["0B"] = [
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 0.64, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 1.5 ],
    [ 0.25, 0.3, 0.41, 0.52, 1.5 ],
    ]

class Analysis:
    """Class to hold analysis cuts, binning, input files, and trigger info"""
    def __init__(self, region, tag, njetsMin=-1, njetsMax=-1, nbMin=-1, nbMax=-1):
        self.region = region
        self.tag = tag
        self.njetsMin = njetsMin
        self.njetsMax = njetsMax
        self.nbMin = nbMin
        self.nbMax = nbMax

        if tag == "Razor2015":
            self.lumi = 2300
        elif tag == "Razor2016":
            #self.lumi = 589
            #self.lumi = 804
            self.lumi = 2600
        else:
            sys.exit("Error: tag"+tag+"is not supported!")

        #load information for analysis
        self.getAnalysis()
        #weight options
        self.getWeightOpts()

    def getAnalysis(self):
        btag = str( max(0, self.nbMin) )+"B"
        if self.region == "SingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["TTJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = razorBinning["SingleLepton"]
            self.unrollBins = (xbinsSignal["TTJetsSingleLepton"]["0B"], colsSignal["TTJetsSingleLepton"]["0B"])

        elif self.region == "SingleLeptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoW"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLeptonInv"][self.tag]
            self.samples = razorSamples["WJetsSingleLeptonInv"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLeptonInv"]
            self.cuts = razorCuts["SingleLeptonInv"]
            self.binning = razorBinning["SingleLeptonInv"]
            self.unrollBins = (xbinsSignal["WJetsSingleLeptonInv"]["0B"], colsSignal["WJetsSingleLeptonInv"]["0B"])

        elif self.region == "TTJetsSingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["TTJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["TTJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = razorBinning["SingleLepton"]
            self.unrollBins = (xbinsSignal["TTJetsSingleLepton"]["0B"], colsSignal["TTJetsSingleLepton"]["0B"])

        elif self.region == "WJetsSingleLepton":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLepton"][self.tag]
            self.samples = razorSamples["WJetsSingleLepton"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLepton"]
            self.cuts = razorCuts["SingleLepton"]
            self.binning = razorBinning["SingleLepton"]
            self.unrollBins = (xbinsSignal["WJetsSingleLepton"]["0B"], colsSignal["WJetsSingleLepton"]["0B"])

        elif self.region == "WJetsSingleLeptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoW"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["SingleLeptonInv"][self.tag]
            self.samples = razorSamples["WJetsSingleLeptonInv"]
            self.samplesReduced = razorSamplesReduced["WJetsSingleLeptonInv"]
            self.cuts = razorCuts["SingleLeptonInv"]
            self.binning = razorBinning["SingleLeptonInv"]
            self.unrollBins = (xbinsSignal["WJetsSingleLeptonInv"]["0B"], colsSignal["WJetsSingleLeptonInv"]["0B"])

        elif self.region == "DYJetsDileptonInv":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoZ"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["DileptonInv"][self.tag]
            self.samples = razorSamples["DYJetsDileptonInv"]
            self.samplesReduced = razorSamplesReduced["DYJetsDileptonInv"]
            self.cuts = razorCuts["DYJetsDileptonInv"]
            self.binning = razorBinning["DYJetsDileptonInv"]
            self.unrollBins = (None,None)

        elif self.region == "DYJetsDileptonInvMultiJet":
            self.trigType = "SingleLepton"
            self.jetVar = "NJets_NoZ"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["DileptonInv"][self.tag]
            self.samples = razorSamples["DYJetsDileptonInv"]
            self.samplesReduced = razorSamplesReduced["DYJetsDileptonInv"]
            self.cuts = razorCuts["DYJetsDileptonInv"]
            self.binning = razorBinning["DYJetsDileptonInvReduced"]
            self.unrollBins = (None,None)

        elif self.region == "VetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["VetoLepton"][self.tag]
            self.samples = razorSamples["VetoLepton"]
            self.samplesReduced = razorSamplesReduced["VetoLepton"]
            self.cuts = razorCuts["VetoLepton"]
            self.binning = razorBinning["VetoLepton"]
            self.unrollBins = (None,None)

        elif self.region == "MultiJetForVetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["MultiJetForVetoLepton"]
            self.binning = razorBinning["SignalRegionForVetoLepton"]
            self.unrollBins = (None,None)
        
        elif self.region == "DiJetForVetoLeptonControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["DiJetForVetoLepton"]
            self.binning = razorBinning["SignalRegionForVetoLepton"]
            self.unrollBins = (None,None)

        elif self.region == "VetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "NJets40"
            self.bjetVar = "NBJetsMedium"
            self.filenames = razorNtuples["VetoTau"][self.tag]
            self.samples = razorSamples["VetoTau"]
            self.samplesReduced = razorSamplesReduced["VetoTau"]
            self.cuts = razorCuts["VetoTau"]
            self.binning = razorBinning["VetoTau"]
            self.unrollBins = (None,None)

        elif self.region == "MultiJetForVetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["MultiJetForVetoTau"]
            self.binning = razorBinning["SignalRegionForVetoTau"]
            self.unrollBins = (None,None)

        elif self.region == "DiJetForVetoTauControlRegion":
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts["DiJetForVetoTau"]
            self.binning = razorBinning["SignalRegionForVetoTau"]
            self.unrollBins = (None,None)

        elif self.region in electronRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalElectron"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning["SignalLeptonic"]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in muonRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalMuon"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning["SignalLeptonic"]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in leptonRazorBoxes:
            self.trigType = "SingleLepton"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalLepton"][self.tag]
            self.samples = razorSamples["SignalLeptonic"]
            self.samplesReduced = razorSamplesReduced["SignalLeptonic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning["SignalLeptonic"]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        elif self.region in hadronicRazorBoxes:
            self.trigType = "Razor"
            self.jetVar = "nSelectedJets"
            self.bjetVar = "nBTaggedJets"
            self.filenames = razorNtuples["SignalHadronic"][self.tag]
            self.samples = razorSamples["SignalHadronic"]
            self.samplesReduced = razorSamplesReduced["SignalHadronic"]
            self.cuts = razorCuts[self.region]
            self.binning = razorBinning["SignalHadronic"]
            self.unrollBins = (xbinsSignal[self.region][btag], colsSignal[self.region][btag])

        else: 
            sys.exit("Error: analysis region "+self.region+" is not implemented!")

        #add jet and bjet and trigger requirements to cuts
        if self.region in razorBoxes: self.cuts = appendNoiseFilters(self.cuts)
        self.addJetCuts()
        self.addTriggerCuts()

    def addJetCuts(self):
        if self.njetsMin >= 0:
            self.cuts = self.cuts + ' && %s >= %d'%(self.jetVar,self.njetsMin)
        if self.njetsMax >= 0:
            self.cuts = self.cuts + ' && %s <= %d'%(self.jetVar,self.njetsMax)
        if self.nbMin >= 0:
            self.cuts = self.cuts + ' && %s >= %d'%(self.bjetVar,self.nbMin)
        if self.nbMax >= 0:
            self.cuts = self.cuts + ' && %s <= %d'%(self.bjetVar,self.nbMax)

    def addTriggerCuts(self):
        self.trigInfoMC = TriggerUtils(self.trigType, self.tag, isData=False)
        self.trigInfoData = TriggerUtils(self.trigType, self.tag, isData=True)
        self.cutsMC = self.trigInfoMC.appendTriggerCuts( self.cuts )+' && TMath::Finite(weight) && (!TMath::IsNaN(weight))'
        self.cutsData = self.trigInfoData.appendTriggerCuts( self.cuts )

    def getWeightOpts(self):
        self.weightOpts = razorWeightOpts[self.tag]
        self.weightFiles = {}
        self.weightHists = {}
        for wType,wInfo in razorWeightHists[self.tag].iteritems():
            self.weightFiles[wType] = rt.TFile.Open(wInfo[0])
            self.weightHists[wType] = self.weightFiles[wType].Get(wInfo[1])
