import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
import macro.macro as macro
from macro.razorAnalysis import razorCuts
from macro.razorWeights import loadScaleFactorHists, invertHistogram
from macro.razorMacros import runFitAndToys, makeControlSampleHists

LUMI = 2185 #in /pb
MCLUMI = 1 

SAMPLES_HADRONIC = ["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets2L", "TTJets1L"]
SAMPLES_LEPTONIC = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"]
SAMPLES = { "MultiJet":SAMPLES_HADRONIC, "MuMultiJet":SAMPLES_LEPTONIC, "EleMultiJet":SAMPLES_LEPTONIC }
BOXES = ["MultiJet", "MuMultiJet", "EleMultiJet"]

DIR_MC = "Backgrounds"
DIR_DATA = "Backgrounds"
#DIR_MC = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108/"
#DIR_DATA = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106"
DATA_NAMES={
    'MultiJet':'RazorInclusive_HTMHT_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered',
    'EleMultiJet':'RazorInclusive_SingleElectron_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered',
    'MuMultiJet':'RazorInclusive_SingleMuon_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered',
    }
FILENAMES_MC = {
        "TTJets1L"    : DIR_MC+"/"+"FullRazorInclusive_TTJets1L_1pb_weighted.root",
        "TTJets2L"    : DIR_MC+"/"+"FullRazorInclusive_TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_1pb_weighted.root",
        "WJets"     : DIR_MC+"/"+"FullRazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root",
        "SingleTop" : DIR_MC+"/"+"FullRazorInclusive_SingleTop_1pb_weighted.root",
        "Other" : DIR_MC+"/"+"FullRazorInclusive_Other_1pb_weighted.root",
        "DYJets"     : DIR_MC+"/"+"FullRazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted.root",
        "ZInv"     : DIR_MC+"/"+"FullRazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root",
        "QCD"       : DIR_DATA+'/'+DATA_NAMES["MultiJet"]+'.root' #data-driven QCD prediction for MultiJet
        }
FILENAMES = {name:copy.copy(FILENAMES_MC) for name in BOXES}
for name in BOXES: FILENAMES[name]["Data"] = DIR_DATA+'/'+DATA_NAMES[name]+'.root'

config = "config/run2_20151108_Preapproval_2b3b_data.config"
FIT_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FitResults/ResultForDecemberJamboree2015/Data_2185ipb"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_MultiJet.root",
        "MuMultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_EleMultiJet.root",
        }

weightOpts = []
commonShapeErrors = [('singletopnorm',"SingleTop"),('othernorm',"Other"),('qcdnorm','QCD'),'btag','pileup','bmistag','facscale','renscale','facrenscale']
commonShapeErrors += [('btagcrosscheckrsq',['TTJets1L','TTJets2L','WJets']),('btagcrosscheckmr',['TTJets1L','TTJets2L','WJets']),('sfsyszinv',['ZInv']),('zllcrosscheck',['ZInv']),'jes','ees','mes',('ttcrosscheck',['TTJets2L']),('sfsysttjets',['TTJets1L','TTJets2L']),('sfsyswjets',['WJets'])]
lepShapeErrors = commonShapeErrors+['tightmuoneff','tighteleeff','muontrig','eletrig']
hadShapeErrors = commonShapeErrors+['sfsysvetolep','sfsysvetotau','mteff','dphieff','vetomuoneff','vetoeleeff']
shapes = { 'MultiJet':hadShapeErrors, 'MuMultiJet':lepShapeErrors, 'EleMultiJet':lepShapeErrors }

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

#scale factor file names
sfdir = "data/ScaleFactors/RazorMADD2015/"
sfFile = sfdir+'/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root'
vetolepFile = sfdir+'/RazorVetoLeptonCrossCheck.root'
vetotauFile = sfdir+'/RazorVetoTauCrossCheck.root'
ttFile = sfdir+'/TTBarDileptonSystematic.root'
dyFile = sfdir+'/RazorDYJetsDileptonInvCrossCheck.root'
btagFile = sfdir+'/RazorBTagClosureTests.root'

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
    parser.add_argument('--no-sys', help="no shape unncertainties or cross check systematics", action="store_true", dest='noSys')
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
    sfNames={
            "ZInv":"GJetsInv",
            "TTJets1L":"TTJets",
            "TTJets2L":"TTJets",
            "ZInvUp":"WJetsInv" #interpolate between GJets and WJets estimates for ZInv scale factors
            }
    sfHists = loadScaleFactorHists(sfFilename=sfFile, processNames=SAMPLES_HADRONIC+['ZInvUp'], scaleFactorNames=sfNames, debugLevel=debugLevel)
    for name in sfHists: assert sfHists[name]

    #get 'down' histogram for wjetsinv/gjets scale factor comparison
    sfHists['ZInvDown'] = sfHists['ZInv'].Clone('WJetsInvScaleFactorsDownFromGJetsInv')
    #TODO: interpolate down from GJets scale factor hist using WJetsInv scale factors (right now hists have diff. sizes)

    #get veto lepton and tau scale factor histograms
    vnames = ['', 'Up', 'Down', 'MTUp', 'MTDown', 'DPhiUp', 'DPhiDown']
    vlfile = rt.TFile.Open(vetolepFile)
    assert vlfile
    vtfile = rt.TFile.Open(vetotauFile)
    assert vtfile
    vlhists = { 'VetoLepton'+n:vlfile.Get('VetoLeptonScaleFactors'+n) for n in vnames }
    vthists = { 'VetoTau'+n:vtfile.Get('VetoTauScaleFactors'+n) for n in vnames }
    for n in vnames: 
        assert vlhists['VetoLepton'+n]
        assert vthists['VetoTau'+n]
    sfHists.update(vlhists)
    sfHists.update(vthists)
    #get DYJets and TTBar Dilepton cross check scale factor histograms
    ttTFile = rt.TFile.Open(ttFile)
    assert ttTFile
    sfHists['TTJetsDileptonUp'] = ttTFile.Get('TTBarDileptonSystematic')
    assert sfHists['TTJetsDileptonUp']
    #convert to correct SF histogram format
    for nb in range(sfHists['TTJetsDileptonUp'].GetSize()+1):
        sfHists['TTJetsDileptonUp'].SetBinContent( nb, sfHists['TTJetsDileptonUp'].GetBinContent(nb)+1.0 )
    #get 'down' version of histogram
    sfHists['TTJetsDileptonDown'] = invertHistogram(sfHists['TTJetsDileptonUp'])
    dyTFile = rt.TFile.Open(dyFile)
    assert dyTFile
    sfHists['DYJetsInvUp'] = dyTFile.Get('DYJetsDileptonInvCrossCheckScaleFactors')
    assert sfHists['DYJetsInvUp']
    sfHists['DYJetsInvDown'] = invertHistogram(sfHists['DYJetsInvUp'])
    btagTFile = rt.TFile.Open(btagFile)
    for b in range(4):
        bs = str(b)
        sfHists['MR'+bs+'BUp'] = btagTFile.Get('OneLepton'+bs+'BMRScaleFactors')
        sfHists['Rsq'+bs+'BUp'] = btagTFile.Get('OneLepton'+bs+'BRsqScaleFactors')
        assert sfHists['MR'+bs+'BUp']
        assert sfHists['Rsq'+bs+'BUp']
        sfHists['MR'+bs+'BDown'] = invertHistogram(sfHists['MR'+bs+'BUp'])
        sfHists['Rsq'+bs+'BDown'] = invertHistogram(sfHists['Rsq'+bs+'BUp'])

    auxSFs = { 
        "VetoLepton":("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13"), 
        "VetoTau":("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 15")
        }

    #estimate yields in signal region
    for boxName in BOXES:

        #apply options
        blindBinsToUse = blindBins[boxName]
        if args.unblind: blindBinsToUse = None
        samplesToUse = SAMPLES[boxName]
        if args.noMC: samplesToUse = []
        if samplesToUse is None or len(samplesToUse) == 0:
            filesToUse = {"Data":FILENAMES[boxName]["Data"]}
        else:
            filesToUse = FILENAMES[boxName]
        if args.noData: 
            del filesToUse['Data']
        shapesToUse = shapes
        if args.noSys:
            shapesToUse = { "MultiJet":[], "MuMultiJet":[], "EleMultiJet":[] } 

        #apply veto lepton correction only to MultiJet box
        if boxName == "MultiJet":
            auxSFsToUse = auxSFs
        else:
            auxSFsToUse = {}

        #loop over btag bins
        #btaglist = [0]
        btaglist = [0,1,2,3]
        for btags in btaglist:
            print "\n---",boxName,"Box,",btags,"B-tags ---"
            #get correct b-tag closure test histogram
            sfHists['MRBUp'] = sfHists['MR'+str(btags)+'BUp']
            sfHists['MRBDown'] = sfHists['MR'+str(btags)+'BDown']
            sfHists['RsqBUp'] = sfHists['Rsq'+str(btags)+'BUp']
            sfHists['RsqBDown'] = sfHists['Rsq'+str(btags)+'BDown']
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
            #if not args.noFit and not os.path.isfile(TOYS_FILES[boxName]):
            #    print "Fit file",TOYS_FILES[boxName],"not found, trying to recreate it"
            #    runFitAndToys(FIT_DIR, boxName, LUMI, DATA_NAMES[boxName], DIR_DATA, config=config, sideband=do#Sideband)
            #    #check
            #    if not os.path.isfile(TOYS_FILES[boxName]):
            #        print "Error creating fit file",TOYS_FILES[boxName]
            #        sys.exit()
            makeControlSampleHists(extboxName, 
                    filenames=filesToUse, samples=samplesToUse, 
                    cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                    bins=binning[boxName], lumiMC=MCLUMI, lumiData=LUMI, 
                    weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                    weightOpts=weightOpts, shapeErrors=shapesToUse[boxName], 
                    fitToyFiles=TOYS_FILES, boxName=boxName, blindBins=blindBinsToUse,
                    btags=nBtags, debugLevel=debugLevel, auxSFs=auxSFsToUse, dataDrivenQCD=True, printdir=dirName, plotOpts=plotOpts)
