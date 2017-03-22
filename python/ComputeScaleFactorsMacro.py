#!/usr/bin/env python
import sys, os, argparse
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--tag", help="Analysis tag, e.g. Razor2015", default="Razor2016")
    parser.add_argument("--no-save", dest="noSave", action="store_true", help="Do not save SFs or histograms")
    parser.add_argument('--no-fill', help="dry run -- do not fill histograms", action="store_true", 
            dest='noFill')
    parser.add_argument('--delta-phi-cut', help='cut on delta phi variable', action='store_true',
            dest='deltaPhiCut')
    parser.add_argument('--njets80-cut', help='require two jets of 80 GeV', action='store_true',
            dest='njets80Cut')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag

    #initialize
    sfHists = {}
    plotOpts = { 'comment':False, "SUS15004CR":True }

    regionsOrder = ["TTJetsSingleLepton", "WJetsSingleLepton", "WJetsSingleLeptonInv"]
    if tag != "Razor2015":
        regionsOrder.insert(0, "GJetsInv")
    regions = {
            "TTJetsSingleLepton":Analysis("TTJetsSingleLepton",tag=tag,nbMin=1),
            "WJetsSingleLepton":Analysis("WJetsSingleLepton",tag=tag,nbMax=0),
            "WJetsSingleLeptonInv":Analysis("WJetsSingleLeptonInv",tag=tag,nbMax=0),
            "GJetsInv":Analysis("GJetsInv",tag=tag)
            }

    for region in regionsOrder:
        analysis = regions[region]
        process = region.replace('SingleLepton','')
        (xbins,cols) = analysis.unrollBins
        #make output directory
        outdir = 'Plots/'+tag+'/'+region
        if args.noSave:
            outdir += '_Test'
        if args.deltaPhiCut:
            outdir += '_DPhiCut'
        if args.njets80Cut:
            outdir += '_NJets80Cut'
        os.system('mkdir -p '+outdir)
        #get correct variable names
        sfVars = ("MR","Rsq")
        if region == "GJetsInv": sfVars = ("MR_NoPho","Rsq_NoPho")
        elif "Inv" in region: sfVars = ("MR_NoW", "Rsq_NoW")
        #QCD estimate for photon region
        if region == 'GJetsInv':
            dataDrivenQCD = True
            auxSFs = razorWeights.getPhotonPuritySFs()
            razorWeights.loadPhotonPurityHists(sfHists, tag, debugLevel)
        else:
            dataDrivenQCD = False
            auxSFs = {}
        #optionally cut on dPhi
        if args.deltaPhiCut:
            dPhiVar = "dPhiRazor"
            if region == 'GJetsInv':
                dPhiVar += '_NoPho'
            analysis.cutsData += " && abs(%s) < 2.8"%(dPhiVar)
            analysis.cutsMC += " && abs(%s) < 2.8"%(dPhiVar)
        if args.njets80Cut:
            njets80Var = "NJets80"
            if region == 'GJetsInv':
                njets80Var += '_NoPho'
            analysis.cutsData += " && %s >= 2"%(njets80Var)
            analysis.cutsMC += " && %s >= 2"%(njets80Var)
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, 
                sfHists=sfHists, sfVars=sfVars, printdir=outdir, debugLevel=debugLevel, 
                auxSFs=auxSFs, noFill=args.noFill, dataDrivenQCD=dataDrivenQCD) 
        #compute scale factors
        appendScaleFactors(process, hists, sfHists, lumiData=analysis.lumi, th2PolyXBins=xbins, 
                th2PolyCols=cols, debugLevel=debugLevel, var=sfVars, printdir=outdir)
        #export histograms
        if not args.noSave:
            macro.exportHists(hists, outFileName='controlHistograms'+region+'.root', 
                outDir=outdir, debugLevel=debugLevel)

    #write scale factors
    if not args.noSave:
        outname = "data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag)
        if args.deltaPhiCut:
            outname = outname.replace('.root','_DPhiCut.root')
        if args.njets80Cut:
            outname = outname.replace('.root','_NJets80Cut.root')
        outfile = rt.TFile(outname, "RECREATE")
        for name in sfHists:
            print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
            sfHists[name].Write(sfHists[name].GetName().replace("Poly",""))
        outfile.Close()
