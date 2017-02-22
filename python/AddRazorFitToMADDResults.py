import sys
import argparse
import ROOT as rt

import macro.macro as macro
from macro.razorAnalysis import razorBinning, razorFitFiles
from macro.razorMacros import import2DRazorFitHistograms

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", default="Razor2016_MoriondRereco",
            help="Analysis tag, e.g. Razor2015")
    parser.add_argument('--box', help="choose a box", required=True)
    args = parser.parse_args()
    
    c = rt.TCanvas("c","c",400,300) # dummy canvas
    bins = razorBinning[args.box]
    btagsMax = 3
    if args.box in ['DiJet','LeptonJet']:
        btagsMax = 2
    for btags in range(btagsMax+1):
        inDir = "Plots/%s/%s%dB"%(args.tag, args.box, btags)
        inFile = "razorHistograms%s%dB.root"%(args.box, btags)
        hists = macro.importHists(inDir+'/'+inFile)
        if 'Fit' in hists:
            sys.exit("Fit histograms already exist in file!")
        import2DRazorFitHistograms(hists, bins, 
                razorFitFiles[args.tag][args.box], c, 
                "Data", btags, btagsMax, noStat=True)
        macro.exportHists(hists, outFileName=inFile, outDir=inDir)
