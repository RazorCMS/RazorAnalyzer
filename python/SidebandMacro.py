import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
from macro.razorAnalysis import razorCuts, xbinsSignal, colsSignal
from macro.razorWeights import loadScaleFactorHists, invertHistogram
from macro.razorMacros import runFitAndToys, makeControlSampleHists
import macro.macro as macro

LUMI = 2185 #in /pb
MCLUMI = 1 

SAMPLES_HADRONIC = ["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets2L", "TTJets1L"]
SAMPLES_LEPTONIC = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"]
SAMPLES = { "MultiJet":SAMPLES_HADRONIC, "MuMultiJet":SAMPLES_LEPTONIC, "EleMultiJet":SAMPLES_LEPTONIC }
BOXES = ["MultiJet", "MuMultiJet", "EleMultiJet"]

DIR_MC = "Backgrounds"
DIR_DATA = "Backgrounds"
#DIR_MC = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p23_Background_20160108/"
#DIR_DATA = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106" #old data 
#DIR_DATA = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForMoriond20160119/RazorSkim"
DATA_NAMES={
    'MultiJet':'RazorInclusive_HTMHT_Run2015D_GoodLumiGolden_RazorSkim_Filtered',
    'EleMultiJet':'RazorInclusive_SingleElectron_Run2015D_GoodLumiGolden_RazorSkim_Filtered',
    'MuMultiJet':'RazorInclusive_SingleMuon_Run2015D_GoodLumiGolden_RazorSkim_Filtered',
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
#FIT_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FitResults/ResultForDecemberJamboree2015/Data_2185ipb" #old fits
FIT_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FitResults/ResultForMoriond2016"
FULL_FIT_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FitResults/ResultForMoriond2016/Full"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_MultiJet.root",
        "MuMultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_EleMultiJet.root",
        }
FULL_TOYS_FILES = {
        #"MultiJet":FULL_FIT_DIR+"/toys_Bayes_varyN_noStat_MultiJet.root",
        #"MuMultiJet":FULL_FIT_DIR+"/toys_Bayes_varyN_noStat_MuMultiJet.root",
        "EleMultiJet":FULL_FIT_DIR+"/toys_Bayes_noStat_EleMultiJet.root",
        }

weightOpts = []
commonShapeErrors = [('singletopnorm',"SingleTop"),('othernorm',"Other"),('qcdnorm','QCD'),'btag','pileup','bmistag','facscale','renscale','facrenscale']
commonShapeErrors += [('btaginvcrosscheck',['ZInv']),('btagcrosscheckrsq',['TTJets1L','TTJets2L','WJets']),('btagcrosscheckmr',['TTJets1L','TTJets2L','WJets']),('sfsyszinv',['ZInv']),('zllcrosscheck',['ZInv']),'jes','ees','mes',('ttcrosscheck',['TTJets2L']),('sfsysttjets',['TTJets1L','TTJets2L']),('sfsyswjets',['WJets'])]
lepShapeErrors = commonShapeErrors+['tightmuoneff','tighteleeff','muontrig','eletrig']
hadShapeErrors = commonShapeErrors+['sfsysvetoleppt','sfsysvetotaupt','mteffpt','dphieffpt','vetolepetacrosscheck','vetotauetacrosscheck','vetomuoneff','vetoeleeff']
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
gjetsupdownFile = sfdir+'/RazorScaleFactors_Inclusive_CorrectedToMultiJet.root'
vetolepPtFile = sfdir+'/RazorVetoLeptonPtCrossCheck.root'
vetotauPtFile = sfdir+'/RazorVetoTauPtCrossCheck.root'
vetolepEtaFile = sfdir+'/RazorVetoLeptonEtaCrossCheck.root'
vetotauEtaFile = sfdir+'/RazorVetoTauEtaCrossCheck.root'
ttFile = sfdir+'/TTBarDileptonSystematic.root'
dyFile = sfdir+'/RazorDYJetsDileptonInvCrossCheck.root'
btagFile = sfdir+'/RazorBTagClosureTests.root'
invbtagFile = sfdir+'/RazorZNuNuBTagClosureTests.root'

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
    parser.add_argument('--no-qcd', help="do not include QCD prediction", action="store_true", dest='noQCD')
    parser.add_argument('--no-fill', help="dry run -- do not fill histograms", action="store_true", dest='noFill')
    parser.add_argument('--b-inclusive', help="do not bin in number of b-tags", action="store_true", dest='bInclusive')
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--export', action="store_true", help="export histograms instead of making plots")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug

    plotOpts = {"ymin":1e-3}

    doSideband=(not args.full)
    toysToUse = TOYS_FILES
    if not doSideband:
        toysToUse = FULL_TOYS_FILES
        dirName += '_Full'
        plotOpts['sideband'] = False
    else:
        plotOpts['sideband'] = True
    if args.unblind:
        dirName += '_Unblinded'
    if args.noFit: 
        toysToUse = None
        del plotOpts['sideband']
    boxesToUse = BOXES
    if args.box is not None:
        boxesToUse = [args.box]

    #initialize
    weightHists = {}
    sfHists = {}

    #make output directory
    os.system('mkdir -p '+dirName)

    ####LOAD ALL SCALE FACTOR HISTOGRAMS

    #get scale factor histograms
    sfNames={
            "ZInv":"GJetsInv",
            "TTJets1L":"TTJets",
            "TTJets2L":"TTJets",
            }
    sfHists = loadScaleFactorHists(sfFilename=sfFile, processNames=SAMPLES_HADRONIC, scaleFactorNames=sfNames, debugLevel=debugLevel)
    for name in sfHists: assert sfHists[name]
    #get histograms for wjetsinv/gjets scale factor comparison
    gjetsupdownTFile = rt.TFile.Open(gjetsupdownFile)
    sfHists['ZInvUp'] = gjetsupdownTFile.Get('WJetsInvScaleFactors')
    sfHists['ZInvDown'] = gjetsupdownTFile.Get('GJetsInvScaleFactors_Down') #down scale factors are (gjets - (wjets-gjets))
    assert sfHists['ZInvUp']
    assert sfHists['ZInvDown']
    #get veto lepton and tau scale factor histograms
    vnames = ['', 'Up', 'Down', 'MTUp', 'MTDown', 'DPhiUp', 'DPhiDown']
    vlPtFile = rt.TFile.Open(vetolepPtFile)
    assert vlPtFile
    vtPtFile = rt.TFile.Open(vetotauPtFile)
    assert vtPtFile
    vlPtHists = { 'VetoLeptonPt'+n:vlPtFile.Get('VetoLeptonPtScaleFactors'+n) for n in vnames }
    vtPtHists = { 'VetoTauPt'+n:vtPtFile.Get('VetoTauPtScaleFactors'+n) for n in vnames }
    for n in vnames: 
        assert vlPtHists['VetoLeptonPt'+n]
        assert vtPtHists['VetoTauPt'+n]
    sfHists.update(vlPtHists)
    sfHists.update(vtPtHists)
    vlEtaFile = rt.TFile.Open(vetolepEtaFile)
    assert vlEtaFile
    vtEtaFile = rt.TFile.Open(vetotauEtaFile)
    assert vtEtaFile
    sfHists['VetoLeptonEtaUp'] = vlEtaFile.Get('VetoLeptonEtaScaleFactors')
    assert sfHists['VetoLeptonEtaUp']
    sfHists['VetoLeptonEtaDown'] = invertHistogram(sfHists['VetoLeptonEtaUp'])
    sfHists['VetoTauEtaUp'] = vtEtaFile.Get('VetoTauEtaScaleFactors')
    assert sfHists['VetoTauEtaUp']
    sfHists['VetoTauEtaDown'] = invertHistogram(sfHists['VetoTauEtaUp'])
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
    #get ZInv b-tag cross check histogram
    invbtagTFile = rt.TFile.Open(invbtagFile)
    sfHists['ZInvBUp'] = invbtagTFile.Get('ZNuNuBTagClosureSysUnc')
    assert sfHists['ZInvBUp']
    #convert to correct SF histogram format
    for nb in range(sfHists['ZInvBUp'].GetSize()+1):
        sfHists['ZInvBUp'].SetBinContent( nb, sfHists['ZInvBUp'].GetBinContent(nb)+1.0 )
    #get 'down' version of histogram
    sfHists['ZInvBDown'] = invertHistogram(sfHists['ZInvBUp'])

    auxSFs = { 
        "VetoLeptonPt":("leadingGenLeptonPt", "(abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13) && leadingGenLeptonPt > 5"), 
        "VetoTauPt":("leadingGenLeptonPt", "abs(leadingGenLeptonType) == 15 && leadingGenLeptonPt > 20"),
        }

    #estimate yields in signal region
    for boxName in boxesToUse:

        #apply options
        blindBinsToUse = blindBins[boxName]
        if args.unblind: blindBinsToUse = None
        samplesToUse = SAMPLES[boxName]
        if args.noQCD and 'QCD' in samplesToUse:
            samplesToUse.remove('QCD')
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
        if args.bInclusive:
            btaglist = [0]
        else:
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
            unrollBins = (xbinsSignal[boxName][str(btags)+'B'], colsSignal[boxName][str(btags)+'B'])
            hists = makeControlSampleHists(extboxName, 
                    filenames=filesToUse, samples=samplesToUse, 
                    cutsMC=thisBoxCuts, cutsData=thisBoxCuts, 
                    bins=binning[boxName], lumiMC=MCLUMI, lumiData=LUMI, 
                    weightHists=weightHists, sfHists=sfHists, treeName="RazorInclusive", 
                    weightOpts=weightOpts, shapeErrors=shapesToUse[boxName], 
                    fitToyFiles=toysToUse, boxName=boxName, blindBins=blindBinsToUse,
                    btags=nBtags, debugLevel=debugLevel, auxSFs=auxSFsToUse, dataDrivenQCD=True, printdir=dirName, 
                    plotOpts=plotOpts, unrollBins=unrollBins, noFill=args.noFill,
                    makePlots = (not args.export), exportShapeErrs=args.export)

            macro.exportHists(hists, outFileName='razorHistograms'+extboxName+'.root', outDir=dirName, debugLevel=debugLevel)
