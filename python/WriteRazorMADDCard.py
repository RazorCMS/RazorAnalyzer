import sys,os
import argparse
import ROOT as rt

#local imports
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import unrollAndStitch
from RunCombine import exec_me
from DustinTuples2DataCard import uncorrelate, writeDataCard_th1
import macro.macro as macro

LUMI = 2185 #in /pb

SAMPLES_HADRONIC = ["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets2L", "TTJets1L"]
SAMPLES_LEPTONIC = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"]
SAMPLES = { "MultiJet":SAMPLES_HADRONIC, "MuMultiJet":SAMPLES_LEPTONIC, "EleMultiJet":SAMPLES_LEPTONIC }

SIGNAL_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined"
BACKGROUND_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorMADD2015"

config="config/run2_20151108_Preapproval_2b3b_data.config"

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--unblind", help="do not blind signal sensitive region", action='store_true')
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--dir', help="output directory", default="SignalRegionPlots", dest='outDir')
    parser.add_argument('--no-sys',dest="noSys",default=False,action='store_true', help="no systematic templates")
    parser.add_argument('--no-stat',dest="noStat",default=False,action='store_true', help="no MC statistics uncertainty")
    #sigal options
    parser.add_argument('-m','--model', default="T1bbbb", help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=int, help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=int, help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=int, help="mass of LSP")

    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    outDir = args.outDir

    if args.box is None:
        print "Please choose an analysis box with --box"
        sys.exit()
    box = args.box
    boxList = box.split('_') #interpret box1_box2_... as a list of boxes to combine

    #make output directory
    os.system('mkdir -p '+outDir)

    for curBox in boxList:
        #get configuration for this box
        samples = SAMPLES[curBox]
        unrollBins = [(xbinsSignal[curBox][str(btags)+'B'], colsSignal[curBox][str(btags)+'B']) for btags in range(4)]

        #make combined unrolled histograms for background
        backgroundHists = unrollAndStitch(curBox, samples=samples, inDir=BACKGROUND_DIR, outDir=outDir, unrollBins=unrollBins, noSys=args.noSys, addStatUnc=(not args.noStat), debugLevel=debugLevel)

        #treat appropriate uncertainties as uncorrelated bin to bin
        toUncorrelate = ['ttcrosscheck','zllcrosscheck','btagcrosscheckmr','btagcrosscheckrsq','btaginvcrosscheck','vetolepetacrosscheck','vetotauetacrosscheck','vetolepptcrosscheck','vetotauptcrosscheck']
        toUncorrelate += ['stat'+curBox+sample for sample in samples]
        for sys in toUncorrelate:
            if 'stat' in sys:
                if args.noStat: continue
                suppressLevel = 0.1
            else:
                if args.noSys: continue
                suppressLevel = 0.0
            uncorrelate(backgroundHists, sys, suppressLevel=suppressLevel)

        #call SMS template maker
        if 'T1' in args.model:
            modelName = 'SMS-'+args.model+'_'+str(args.mGluino)+'_'+str(args.mLSP)
        else:
            modelName = 'SMS-'+args.model+'_'+str(args.mStop)+'_'+str(args.mLSP)
        signalFilename=SIGNAL_DIR+'/'+modelName+'.root'
        exec_me('python python/SMSTemplates.py --merge-bins -c %s -d %s --lumi %d --box %s %s %s'%(config, outDir, LUMI, curBox, ((args.noSys)*('--no-signal-sys')), signalFilename), False) 
        #load SMS template histograms
        signalHistFilename = '%s/%s_lumi-%.3f_0-3btag_%s.root'%(outDir,modelName,LUMI*1.0/1000,curBox)
        signalHists = macro.importHists(signalHistFilename)
        #update with correct names
        for x,h in signalHists.items():
            h.SetName(h.GetName().replace(curBox+'_'+args.model,modelName))
            signalHists[h.GetName()] = signalHists.pop(x)
        #combine signal and background dictionaries
        hists = backgroundHists.copy()
        hists.update(signalHists)

        #write histograms to ROOT file
        outFileName = outDir+'/RazorInclusiveMADD_lumi-%.1f_%s.root'%(LUMI*1.0/1000.,curBox)
        outFile = rt.TFile(outFileName, 'recreate')
        for x,h in hists.iteritems():
            h.Write()
        outFile.Close()
        #write combine card
        cardName = outFileName.replace('.root','.txt')
        writeDataCard_th1(curBox,modelName,cardName,hists,samples)

    #run combine
    if len(boxList) == 1:
        #get card name
        cardName = outDir+'/RazorInclusiveMADD_lumi-%.1f_%s.txt'%(LUMI*1.0/1000.,boxList[0])
        combineName = 'MADD_'+boxList[0]+'_'+modelName
        exec_me('combine -M Asymptotic %s -n %s'%(cardName,combineName), False)
        exec_me('mv higgsCombine%s.Asymptotic.mH120.root %s/'%(combineName,outDir),False)
    elif len(boxList) > 1:
        #get card names
        cardNames = []
        combineName = 'MADD_'+('_'.join(boxList))+'_'+modelName
        for curBox in boxList:
            cardName = 'RazorInclusiveMADD_lumi-%.1f_%s.txt'%(LUMI*1.0/1000.,curBox)
            cardNames.append(cardName)
        #combine cards
        combinedCardName = 'RazorInclusiveMADD_'+('_'.join(boxList))+'.txt'
        exec_me('cd '+outDir+'; combineCards.py '+(' '.join(cardNames))+' > '+combinedCardName+'; cd ..', False)
        #call combine
        exec_me('combine -M Asymptotic '+outDir+'/'+combinedCardName+' -n '+combineName, False)
        #store output in output directory
        exec_me('mv higgsCombine%s.Asymptotic.mH120.root %s/'%(combineName,outDir),False)
