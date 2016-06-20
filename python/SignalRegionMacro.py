import sys,os,argparse,copy
import ROOT as rt

from macro import macro
from macro.razorAnalysis import Analysis
from macro.razorMacros import runFitAndToys, makeControlSampleHists

FIT_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FitResults/ResultForMoriond2016"
FULL_FIT_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FitResults/ResultForMoriond2016/Full"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_MultiJet.root",
        "MuMultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_MuMultiJet.root",
        "EleMultiJet":FIT_DIR+"/toys_Bayes_varyN_noStat_EleMultiJet.root",
        "DiJet":None, "MuJet":None, "EleJet":None, "FourToSixJet":None, "SevenJet":None,
        "LeptonJet":None, "LeptonMultiJet":None,
        }
FULL_TOYS_FILES = {
        "EleMultiJet":FULL_FIT_DIR+"/toys_Bayes_noStat_EleMultiJet.root",
        }

commonShapeErrors = [
        ('singletopnorm',"SingleTop"),
        ('othernorm',"Other"),
        ('qcdnorm','QCD'),
        'btag', 'pileup', 'bmistag', 'facscale', 'renscale', 'facrenscale',
        ('btaginvcrosscheck',['ZInv']),
        ('btagcrosscheckrsq',['TTJets1L','TTJets2L','WJets']),
        ('btagcrosscheckmr',['TTJets1L','TTJets2L','WJets']),
        ('sfstatzinv',['ZInv']),
        ('sfsyszinv',['ZInv']),
        'jes','ees','mes',
        ('ttcrosscheck',['TTJets2L']),
        ('sfstatttjets',['TTJets1L','TTJets2L']),
        ('sfsysttjets',['TTJets1L','TTJets2L']),
        ('sfstatwjets',['WJets']),
        ('sfsyswjets',['WJets'])
        ]
lepShapeErrors = commonShapeErrors+['tightmuoneff','tighteleeff','muontrig','eletrig']
hadShapeErrors = commonShapeErrors+['vetolepptcrosscheck','vetotauptcrosscheck',
        'vetolepetacrosscheck','vetotauetacrosscheck','vetomuoneff','vetoeleeff']
shapes = { 'MultiJet':hadShapeErrors, 'MuMultiJet':lepShapeErrors, 'EleMultiJet':lepShapeErrors, 
           'DiJet':hadShapeErrors, 'MuJet':lepShapeErrors, 'EleJet':lepShapeErrors,
           'FourToSixJet':hadShapeErrors, 'SevenJet':hadShapeErrors }

blindBins = {b:[(x,y) for x in range(2,len(binning[b]["MR"])+1) for y in range(2,len(binning[b]["Rsq"])+1)] for b in binning}

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--unblind", help="do not blind signal sensitive region", action='store_true')
    parser.add_argument('--no-mc', help="do not process MC, do data and fit only", 
            action='store_true', dest="noMC")
    parser.add_argument('--no-fit', help="do not load fit results, process data and MC only", 
            action='store_true', dest='noFit')
    parser.add_argument('--full', help="do full fit (default is sideband)", action='store_true')
    parser.add_argument('--no-data', help="do not process data, do fit and MC only", 
            action='store_true', dest='noData')
    parser.add_argument('--no-sys', help="no shape unncertainties or cross check systematics", 
            action="store_true", dest='noSys')
    parser.add_argument('--no-qcd', help="do not include QCD prediction", action="store_true", 
            dest='noQCD')
    parser.add_argument('--no-fill', help="dry run -- do not fill histograms", action="store_true", 
            dest='noFill')
    parser.add_argument('--no-sfs', help="ignore MC scale factors", action="store_true", 
            dest="noSFs")
    parser.add_argument('--b-inclusive', help="do not bin in number of b-tags", action="store_true",
            dest='bInclusive')
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--btags', type=int, help="choose a number of btags")
    parser.add_argument("--tag", dest="tag", default="Razor2015",
            help="Analysis tag, e.g. Razor2015")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    if tag not in ["Razor2015","Razor2016"]:
        sys.exit("Error: tag "+tag+" not supported!")

    #initialize
    plotOpts = {"SUS15004":True}
    weightHists = {}

    doSideband=(not args.full)
    toysToUse = TOYS_FILES
    dirSuffix = ""
    if not doSideband:
        toysToUse = FULL_TOYS_FILES
        dirSuffix += '_Full'
        plotOpts['sideband'] = False
    else:
        plotOpts['sideband'] = True
    if args.unblind:
        dirSuffix += '_Unblinded'
    if args.noFit: 
        toysToUse = {}
        del plotOpts['sideband']
    boxesToUse = ["MultiJet", "MuMultiJet", "EleMultiJet", "DiJet", "MuJet", "EleJet"]
    if args.box is not None:
        boxesToUse = [args.box]
    if args.btags is not None:
        btaglist = [args.btags]
    else:
        btaglist = [0,1,2,3]
    if args.noSFs:
        dirSuffix += '_NoSFs'

    regionsOrder = []
    regions = {}
    for box in boxesToUse:
        for btags in btaglist:
            #deal with last b-tag bin
            if box in ["DiJet", "MuJet", "EleJet", "LeptonJet"]:
                if btags == 3: continue
                elif btags >= 2:
                    nbMax = -1 #no upper limit
                else:
                    nbMax = btags #exclusive
            else:
                if btags >= 3:
                    nbMax = -1 #no upper limit
                else:
                    nbMax = btags #exclusive
            #define analysis region
            extBox = box+str(btags)+"B"
            regionsOrder.append(extBox)
            regions[extBox] = Analysis(box, tag=tag, nbMin=btags, nbMax=nbMax)

    ####LOAD ALL SCALE FACTOR HISTOGRAMS

    #scale factor file names
    sfdir = "data/ScaleFactors/RazorMADD2015/"
    sfFile = sfdir+'/RazorScaleFactors_%s.root'%(tag)
    sfFile_nJets = sfdir+'/RazorNJetsScaleFactors_%s.root'%(tag)
    vetolepFile = sfdir+'/RazorVetoLeptonCrossCheck_%s.root'%(tag)
    ttFile = sfdir+'/TTBarDileptonSystematic_%s.root'%(tag)
    dyFile = sfdir+'/RazorDYJetsDileptonInvCrossCheck_%s.root'%(tag)
    btagFile = sfdir+'/RazorBTagClosureTests_%s.root'%(tag)

    #get MR-Rsq scale factor histograms
    sfNames={
            "ZInv":"GJetsInv",
            "TTJets1L":"TTJets",
            "TTJets2L":"TTJets",
            }
    processNames = regions.itervalues.next().samples
    sfHists = loadScaleFactorHists(sfFilename=sfFile, processNames=processNames, 
            scaleFactorNames=sfNames, debugLevel=debugLevel)
    #reopen the file and grab the ZNuNu up/down histograms
    #down scale factors are (gjets - (wjets-gjets))
    sfTFile = rt.TFile.Open(gjetsupdownFile)
    sfHists['ZInvUp'] = sfTFile.Get('WJetsInvScaleFactors')
    sfHists['ZInvDown'] = sfTFile.Get('GJetsInvScaleFactors_Down') 
    #get njets scale factor histogram
    sfNJetsFile = rt.TFile.Open(sfFile_nJets)
    sfHists['NJets'] = sfNJetsFile.Get("NJetsCorrectionScaleFactors")
    #get veto lepton and tau scale factor histograms
    vlFile = rt.TFile.Open(vetolepFile)
    sfHists['VetoLeptonPtUp'] = vlFile.Get('VetoLeptonPtScaleFactors')
    sfHists['VetoLeptonPtDown'] = invertHistogram(sfHists['VetoLeptonPtUp'])
    sfHists['VetoTauPtUp'] = vlFile.Get('VetoTauPtScaleFactors')
    sfHists['VetoTauPtDown'] = invertHistogram(sfHists['VetoTauPtUp'])
    sfHists['VetoLeptonEtaUp'] = vlFile.Get('VetoLeptonEtaScaleFactors')
    sfHists['VetoLeptonEtaDown'] = invertHistogram(sfHists['VetoLeptonEtaUp'])
    sfHists['VetoTauEtaUp'] = vlFile.Get('VetoTauEtaScaleFactors')
    sfHists['VetoTauEtaDown'] = invertHistogram(sfHists['VetoTauEtaUp'])
    #get DYJets and TTBar Dilepton cross check scale factor histograms
    ttTFile = rt.TFile.Open(ttFile)
    sfHists['TTJetsDileptonUp'] = ttTFile.Get('TTBarDileptonSystematic')
    #convert to correct SF histogram format
    for nb in range(sfHists['TTJetsDileptonUp'].GetSize()+1):
        sfHists['TTJetsDileptonUp'].SetBinContent( nb, 
                sfHists['TTJetsDileptonUp'].GetBinContent(nb)+1.0 )
    #get 'down' version of histogram
    sfHists['TTJetsDileptonDown'] = invertHistogram(sfHists['TTJetsDileptonUp'])
    dyTFile = rt.TFile.Open(dyFile)
    sfHists['DYJetsInvUp'] = dyTFile.Get('DYJetsDileptonInvCrossCheckScaleFactors')
    sfHists['DYJetsInvDown'] = invertHistogram(sfHists['DYJetsInvUp'])
    btagTFile = rt.TFile.Open(btagFile)
    for b in range(4):
        bs = str(b)
        sfHists['MR'+bs+'BUp'] = btagTFile.Get('OneLepton'+bs+'BMRScaleFactors')
        sfHists['Rsq'+bs+'BUp'] = btagTFile.Get('OneLepton'+bs+'BRsqScaleFactors')
        sfHists['MR'+bs+'BDown'] = invertHistogram(sfHists['MR'+bs+'BUp'])
        sfHists['Rsq'+bs+'BDown'] = invertHistogram(sfHists['Rsq'+bs+'BUp'])
    #get ZInv b-tag cross check histogram
    invbtagTFile = rt.TFile.Open(invbtagFile)
    sfHists['ZInvBUp'] = invbtagTFile.Get('ZNuNuBTagClosureSysUnc')
    #convert to correct SF histogram format
    for nb in range(sfHists['ZInvBUp'].GetSize()+1):
        sfHists['ZInvBUp'].SetBinContent( nb, sfHists['ZInvBUp'].GetBinContent(nb)+1.0 )
    #get 'down' version of histogram
    sfHists['ZInvBDown'] = invertHistogram(sfHists['ZInvBUp'])

    auxSFs = {"NJets":("nSelectedJets","1")} #do not correct veto lepton pt or eta

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
        shapesToUse = copy.copy(shapes[boxName])
        if args.noSys:
            shapesToUse = []
        if (boxName not in toysToUse or toysToUse[boxName] is None) and 'sideband' in plotOpts:
            del plotOpts['sideband']

        sfHistsToUse = sfHists
        auxSFsToUse = auxSFs
        #if boxName == 'SevenJet':
            #sfHistsToUse = sfHists7Jet

        #disable scale factors option
        if args.noSFs:
            print "Ignoring all scale factor histograms and uncertainties from scale factor cross checks."
            sfHistsToUse = {}
            auxSFsToUse = {}
            toRemove = ['btaginvcrosscheck','btagcrosscheckrsq','btagcrosscheckmr','sfsyszinv','ttcrosscheck','zllcrosscheck','sfsysttjets','sfsyswjets','vetolepptcrosscheck','vetotauptcrosscheck','vetolepetacrosscheck','vetotauetacrosscheck']
            #remove scale factor cross check uncertainties
            shapesToUse = [s for s in shapesToUse if s not in toRemove]
            #this removes scale factor uncertainties that are listed as tuples
            shapesToUse = [s for s in shapesToUse if not (hasattr(s, '__getitem__') and s[0] in toRemove)] 

        #loop over btag bins
        if args.btags is not None:
            btaglist = [args.btags]
        elif args.bInclusive:
            btaglist = [0]
        else:
            btaglist = [0,1,2,3]
        for btags in btaglist:
            print "\n---",boxName,"Box,",btags,"B-tags ---"
            #get correct b-tag closure test histogram
            if not args.noSFs:
                sfHistsToUse['MRBUp'] = sfHistsToUse['MR'+str(btags)+'BUp']
                sfHistsToUse['MRBDown'] = sfHistsToUse['MR'+str(btags)+'BDown']
                sfHistsToUse['RsqBUp'] = sfHistsToUse['Rsq'+str(btags)+'BUp']
                sfHistsToUse['RsqBDown'] = sfHistsToUse['Rsq'+str(btags)+'BDown']
            #get correct cuts string
            thisBoxCuts = razorCuts[boxName]
            if btags >= 3 or args.bInclusive: #inclusive if requested or if we are doing 3B
                thisBoxCuts += " && nBTaggedJets >= "+str(btags)
            else:
                thisBoxCuts += " && nBTaggedJets == "+str(btags)

            if not args.bInclusive:
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
                    weightHists=weightHists, sfHists=sfHistsToUse, treeName="RazorInclusive", 
                    weightOpts=weightOpts, shapeErrors=shapesToUse, 
                    fitToyFiles=toysToUse, boxName=boxName, blindBins=blindBinsToUse,
                    btags=nBtags, debugLevel=debugLevel, auxSFs=auxSFsToUse, dataDrivenQCD=True, printdir=dirName, 
                    plotOpts=plotOpts, unrollBins=unrollBins, noFill=args.noFill,
                    makePlots = (not args.export), exportShapeErrs=args.export, propagateScaleFactorErrs=False)

            macro.exportHists(hists, outFileName='razorHistograms'+extboxName+'.root', outDir=dirName, debugLevel=debugLevel)
