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
MCLUMI = 1 

SAMPLES_HADRONIC = ["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets2L", "TTJets1L"]
SAMPLES_LEPTONIC = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"]
SAMPLES = { "MultiJet":SAMPLES_HADRONIC, "MuMultiJet":SAMPLES_LEPTONIC, "EleMultiJet":SAMPLES_LEPTONIC }

SIGNAL_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined"

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
    parser.add_argument('--dir', help="output directory (should contain the ROOT files with the background histograms)", default="SignalRegionPlots", dest='dirName')
    parser.add_argument('--no-sys',dest="noSys",default=False,action='store_true', help="no systematic templates")
    parser.add_argument('--no-stat',dest="noStat",default=False,action='store_true', help="no MC statistics uncertainty")
    #sigal options
    parser.add_argument('-m','--model', default="T1bbbb", help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=int, help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=int, help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=int, help="mass of LSP")

    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    dirName = args.dirName

    if args.box is None:
        print "Please choose an analysis box with --box"
        sys.exit()
    box = args.box
    #make output directory
    os.system('mkdir -p '+dirName)
    #get configuration for this box
    samples = SAMPLES[box]
    unrollBins = [(xbinsSignal[box][str(btags)+'B'], colsSignal[box][str(btags)+'B']) for btags in range(4)]

    #make combined unrolled histograms for background
    backgroundHists = unrollAndStitch(box, samples=samples, directory=dirName, unrollBins=unrollBins, noSys=args.noSys, addStatUnc=(not args.noStat), debugLevel=debugLevel)

    #treat appropriate uncertainties as uncorrelated bin to bin
    toUncorrelate = ['ttcrosscheck','zllcrosscheck','btagcrosscheckmr','btagcrosscheckrsq','btaginvcrosscheck','vetolepetacrosscheck','vetotauetacrosscheck','vetolepptcrosscheck','vetotauptcrosscheck']
    toUncorrelate += ['stat'+box+sample for sample in samples]
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
    exec_me('python python/SMSTemplates.py --merge-bins -c %s -d %s --lumi %d --box %s %s %s'%(config, dirName, LUMI, box, ((args.noSys)*('--no-signal-sys')), signalFilename), False) 
    #load SMS template histograms
    signalHistFilename = '%s/%s_lumi-%.3f_0-3btag_%s.root'%(dirName,modelName,LUMI*1.0/1000,box)
    signalHists = macro.importHists(signalHistFilename)
    #update with correct names
    for x,h in signalHists.items():
        h.SetName(h.GetName().replace(box+'_'+args.model,modelName))
        signalHists[h.GetName()] = signalHists.pop(x)
    #combine signal and background dictionaries
    hists = backgroundHists.copy()
    hists.update(signalHists)

    #write histograms to ROOT file
    outFileName = dirName+'/RazorInclusiveMADD_lumi-%.1f_%s.root'%(LUMI*1.0/1000.,box)
    outFile = rt.TFile(outFileName, 'recreate')
    for x,h in hists.iteritems():
        h.Write()
    outFile.Close()
    #write combine card
    cardName = outFileName.replace('.root','.txt')
    writeDataCard_th1(box,modelName,cardName,hists,samples)
    #run combine
    exec_me('combine -M Asymptotic '+cardName, False) 
