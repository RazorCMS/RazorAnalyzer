import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from macro import macro
from macro.razorAnalysis import Analysis, razorSignalDirs, signalConfig, controlConfig
from macro.razorMacros import unrollAndStitch, unrollAndStitchFromFiles, getMaxBtags
from RunCombine import exec_me
from SMSTemplates import makeSMSTemplates, signalShapeUncerts, SMSOpts
from DustinTuples2DataCard import uncorrelate, uncorrelateSFs, writeDataCard_th1
from framework import Config
import CheckSignalContamination as contam

BACKGROUND_DIR = "/eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorMADD2016"

def getModelName(model, mass1, mass2):
    return "SMS-%s_%d_%d"%(model, mass1, mass2)

def getBranchingFracsFromModelName(model):
    """Returns (x branching ratio, y branching ratio)"""
    xBR = float(model[model.find('x')+1:
        model.find('y')].replace('p','.'))
    yBR = float(model[model.find('y')+1:].replace(
        'p','.'))
    return xBR, yBR

def getCardName(modelName, box, outDir):
    return outDir+'/RazorInclusiveMADD_%s_%s.txt'%(
                modelName, box)


if __name__ == "__main__":
    rt.gROOT.SetBatch()

    # Verbosity
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true",
            help="display detailed output messages")
    parser.add_argument("-d", "--debug", action="store_true",
            help="display excruciatingly detailed output messages")
    # Basic configuration
    parser.add_argument('--tag', default='Razor2016_MoriondRereco')
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--dir', help="output directory",
            default="SignalRegionPlots", dest='outDir')
    # Customization
    parser.add_argument('--no-limit', dest='noCombine', 
            action='store_true', 
            help='do not call combine, make template histograms only')
    parser.add_argument('--signif', action='store_true', 
            help='compute significance rather than limit')
    parser.add_argument('--no-sys',dest="noSys",default=False,
            action='store_true', help="no systematic templates")
    parser.add_argument('--no-stat',dest="noStat",
            default=False,action='store_true', 
            help="no MC statistics uncertainty")
    parser.add_argument('--no-signal-contam', dest='noSignalContam',
            action='store_true', help='ignore signal contamination')
    parser.add_argument('--fit-sys', dest="addMCVsFit", 
            action='store_true', help="add MC vs fit systematic")
    parser.add_argument('--save-workspace', dest='saveWorkspace', action='store_true',
            help='save combine workspace in output file')
    # Signal model
    parser.add_argument('-m','--model', default="T1bbbb", 
            help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=int, 
            help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=int, 
            help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=int, 
            help="mass of LSP")

    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    outDir = args.outDir
    cfg = Config.Config(signalConfig)

    if args.box is None:
        print "Please choose an analysis box with --box"
        sys.exit()
    box = args.box
    boxList = box.split('_') # list of boxes to combine, separated by _

    # make output directory
    os.system('mkdir -p '+outDir)

    for curBox in boxList:
        btagsMax = getMaxBtags(curBox)

        # retrieve binning and other info
        analyses = []
        unrollBins = []
        for nb in range(btagsMax + 1):
            analyses.append(Analysis(curBox, args.tag, nbMin=nb))
            unrollBins.append(analyses[-1].unrollBins)
        lumi = analyses[0].lumi
        samples = analyses[0].samples
        sfNames = { "TTJets1L":"TTJets", "TTJets2L":"TTJets", 
                "ZInv":"GJetsInv" }
        sfHists = macro.loadScaleFactorHists(
            "data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(args.tag), 
            processNames=["TTJets1L","TTJets2L","WJets","ZInv"], 
            scaleFactorNames=sfNames, debugLevel=debugLevel)

        # assess signal contamination in control regions
        contamHists = None
        if not args.noSignalContam:
            print "Computing signal contamination level in control regions"
            ttContamHist = contam.checkSignalContamination(controlConfig, 
                   outDir=outDir, lumi=lumi, debugLevel=debugLevel,
                   box="TTJetsSingleLeptonForSignal", 
                   model=args.model, mLSP=args.mLSP, 
                   mGluino=args.mGluino, mStop=args.mStop, 
                   mergeBins=True)
            wContamHist = contam.checkSignalContamination(controlConfig,
                   outDir=outDir, lumi=lumi, debugLevel=debugLevel,
                   box="WJetsSingleLeptonForSignal", model=args.model,
                   mLSP=args.mLSP, mGluino=args.mGluino, 
                   mStop=args.mStop, mergeBins=True)
            contamHists = { "TTJets1L":ttContamHist, 
                    "TTJets2L":ttContamHist, "WJets":wContamHist }

        # make combined unrolled histograms for background
        print "Retrieving background histograms from files"
        backgroundHists = unrollAndStitchFromFiles(curBox, 
                samples=samples, inDir=BACKGROUND_DIR, 
                outDir=outDir, unrollBins=unrollBins, noSys=args.noSys, 
                addStatUnc=(not args.noStat), 
                addMCVsFit=args.addMCVsFit, debugLevel=debugLevel)

        if args.noCombine: continue #if not setting a limit, we are done

        # get file name for signal input
        signalDir = razorSignalDirs[args.tag]
        isGluinos = ('T1' in args.model or 'T5' in args.model)
        if isGluinos:
            modelName = getModelName(args.model, args.mGluino, args.mLSP)
        else:
            modelName = getModelName(args.model, args.mStop, args.mLSP)
        signalFilename=signalDir+'/'+modelName+'.root'

        # to modify branching ratios in T1ttbb sample
        xBR = yBR = -1
        if 'T1x' in args.model:
            xBR, yBR = getBranchingFracsFromModelName(args.model)
            signalFilename = signalDir+'/'+getModelName('T1ttbb', 
                    args.mGluino, args.mLSP)+'.root'

        # get correct list of uncertainties for this box
        uncerts = copy.copy(signalShapeUncerts)
        if curBox in ['DiJet','MultiJet']:
            uncertsToRemove = ['tightmuoneff','tighteleeff',
                    'muontrig','eletrig','tightmuonfastsim',
                    'tightelefastsim']
        else:
            uncertsToRemove = ['vetomuoneff','vetoeleeff',
                    'vetomuonfastsim','vetoelefastsim']
        for u in uncertsToRemove:
            if u in uncerts:
                uncerts.remove(u)

        # call SMS template maker
        smsOpts = SMSOpts(xBR=xBR, yBR=yBR, doNPVExtrap=True,
                doGenMetVsPFMet=True)
        signalHists = makeSMSTemplates(curBox, signalFilename,
                uncertainties=uncerts, debugLevel=debugLevel,
                tag=args.tag, opts=smsOpts)
    
        # reduced efficiency method -- corrects for signal contamination
        if not args.noSignalContam:
            beforeHist = signalHists['Signal'].Clone()
            beforeHist.SetLineColor(rt.kRed)
            macro.doDeltaBForReducedEfficiencyMethod(backgroundHists, 
                    signalHists, contamHists, sfHists, 
                    unrollBins=unrollBins, debugLevel=debugLevel)

        # combine signal and background dictionaries
        hists = backgroundHists.copy()
        hists.update(signalHists)

        # treat appropriate uncertainties as uncorrelated bin to bin
        toUncorrelate = ['stat'+curBox+sample for sample in samples]
        for sys in toUncorrelate:
            if 'stat' in sys:
                if args.noStat: continue
                suppressLevel = 0.1
            else:
                if args.noSys: continue
                suppressLevel = 0.0
            uncorrelate(hists, sys, 
                    suppressLevel=suppressLevel)

        # scale factor uncertainties are correlated according to 
        # the bin they are in
        toUncorrelateSF = ['sfstatttjets','sfstatwjets','sfstatzinv']
        for sys in toUncorrelateSF:
            uncorrelateSFs(hists, sys, sfHists, cfg, 
                    curBox, unrollBins=unrollBins)

        # write histograms to ROOT file
        cardName = getCardName(modelName, curBox, outDir)
        outFileName = cardName.replace('.txt', '.root')
        outFile = rt.TFile(outFileName, 'recreate')
        sortedKeys = sorted(hists.keys())
        for key in sortedKeys:
            # b mistag systematic is currently messed up -- ignore
            if 'bmistag' in key.lower():
                del hists[key]
                continue
            hists[key].Write()
        outFile.Close()

        writeDataCard_th1(curBox,cardName,hists,samples)

    if args.noCombine: 
        sys.exit()

    # get combine parameters
    if args.signif:
        combineMethod = 'ProfileLikelihood'
        combineFlags = '--signif -t -1 --toysFreq'
    else:
        combineMethod = 'Asymptotic'
        combineFlags = ''
    if args.saveWorkspace:
        combineFlags += ' --saveWorkspace'

    # run combine 
    if len(boxList) == 1:
        cardName = getCardName(modelName, boxList[0], outDir)
        combineName = 'MADD_'+boxList[0]+'_'+modelName
        exec_me('combine -M %s %s -n %s %s'%(
            combineMethod,cardName,combineName,combineFlags), False)
        exec_me('mv higgsCombine%s.%s.mH120.root %s/'%(
            combineName,combineMethod,outDir),False)
    elif len(boxList) > 1:
        cardNames = []
        combineName = 'MADD_'+('_'.join(boxList))+'_'+modelName
        for curBox in boxList:
            cardName = getCardName(modelName, curBox, outDir)
            cardNames.append(os.path.basename(cardName))
        combinedCardName = ('RazorInclusiveMADD_%s_%s.txt'%(
            modelName, '_'.join(boxList)))
        exec_me('cd '+outDir+'; combineCards.py '+(' '.join(cardNames))
                +' > '+combinedCardName+'; cd ..', False)
        exec_me('combine -M '+combineMethod+' '+outDir+'/'
                +combinedCardName+' -n '+combineName+' '+combineFlags, 
                False)
        exec_me('mv higgsCombine%s.%s.mH120.root %s/'%(
            combineName,combineMethod,outDir),False)
