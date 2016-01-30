import sys,os
import argparse
import ROOT as rt

#local imports
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.razorMacros import unrollAndStitch
from RunCombine import exec_me

LUMI = 2185 #in /pb
MCLUMI = 1 

SAMPLES_HADRONIC = ["Other", "QCD", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets2L", "TTJets1L"]
SAMPLES_LEPTONIC = ["Other", "DYJets", "ZInv", "SingleTop", "WJets", "TTJets1L", "TTJets2L"]
SAMPLES = { "MultiJet":SAMPLES_HADRONIC, "MuMultiJet":SAMPLES_LEPTONIC, "EleMultiJet":SAMPLES_LEPTONIC }

SIGNAL_DIR = {
          'T1bbbb': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForApproval20151208/jobs/combined/',
          'T1tttt': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForApproval20151208/jobs/combined/',
          'T1qqqq': 'root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForApproval20151208/jobs/combined/'
          }

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
    #sigal options
    parser.add_argument('-m','--model', default="T1bbbb", help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=float, help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=float, help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=float, help="mass of LSP")

    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    dirName = args.dirName

    if args.box is None:
        print "Please choose an analysis box with --box"
        sys.exit()

    #make output directory
    os.system('mkdir -p '+dirName)

    #get configuration for this box
    samples = SAMPLES[args.box]
    unrollBins = [(xbinsSignal[args.box][str(btags)+'B'], colsSignal[args.box][str(btags)+'B']) for btags in range(4)]

    #make combined unrolled histograms for background
    unrollAndStitch(args.box, samples=samples, directory=dirName, unrollBins=unrollBins, debugLevel=debugLevel)

    #call SMS template maker
    #if 'T1' in args.model:
    #    signalFilename=SIGNAL_DIR[args.model]+'/SMS-'+args.model+'_'+args.mGluino+'_'+args.mLSP+'.root'
    #else:
    #    signalFilename=SIGNAL_DIR[args.model]+'/SMS-'+args.model+'_'+args.mStop+'_'+args.mLSP+'.root'
    #exec_me('python SMSTemplates.py -c %s -d %s --lumi %d --box %s %s %s'%(config, dirName, LUMI, args.box, ((args.noSignalSys)*('--no-signal-sys')), signalFilename), True) #set to False when ready
