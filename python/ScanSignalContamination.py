import sys
import argparse
import ROOT as rt

#local imports
from framework import Config
from GChiPairs import gchipairs
from SidebandMacro import LUMI
from CheckSignalContamination import checkSignalContamination

def usage():
    print "usage: ScanSignalContamination.py modelName"

if __name__ == '__main__':
    rt.gROOT.SetBatch()

    if len(sys.argv) != 2:
        usage()
        sys.exit()
    model = sys.argv[1]
    pairs = gchipairs(model)

    contam = rt.TH2F("contam", "Signal contamination", 80, 0, 2000, 40, 0, 1000)
    for p in pairs:
        try:
            h = checkSignalContamination("config/run2_20151229_ControlRegion.config", outDir='SMSPlots', lumi=LUMI, box="TTJetsSingleLeptonControlRegion", model=model, mGluino=p[0], mLSP=p[1], mergeBins=True, treeName="RazorInclusive")
            #get maximum signal contamination
            maxContam = 0.0
            for bx in range(1,h.GetSize()):
                maxContam = max(maxContam, h.GetBinContent(bx)) 
            contam.Fill(p[0],p[1],maxContam)
            print p, "signal contamination", maxContam
        except:
            print "Problem getting histogram for",p

    c = rt.TCanvas("c", "c", 800, 600)
    contam.GetXaxis().SetTitle("Gluino mass")
    contam.GetYaxis().SetTitle("LSP mass")
    contam.SetStats(0)
    contam.Draw("colz")
    c.Print("signalContaminationScan"+model+".root")
