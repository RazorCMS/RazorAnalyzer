### Load fit results from toy file and print out the yields bin by bin 

import sys,os
import argparse
import ROOT as rt

#local imports
import macro.macro as macro
from macro.razorMacros import *

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("inFile", help="Toy file used to load fit results")
    parser.add_argument("-c", "--config", help="Config file to use", default="config/run2_sideband.config")
    parser.add_argument("-b", "--box", help="Analysis box", default="MultiJet")
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug

    fitHist = get3DRazorFitHistogram(args.config, args.inFile, args.box, debugLevel=debugLevel)
    makeRazor3DTable(fitHist, args.box)
