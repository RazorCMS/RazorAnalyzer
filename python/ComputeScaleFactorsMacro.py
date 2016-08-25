import sys, os, argparse
import ROOT as rt

from macro import macro
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
    parser.add_argument("--tag", dest="tag", required=True,
                                help="Analysis tag, e.g. Razor2015")
    parser.add_argument("--no-save", dest="noSave", action="store_true", help="Do not save SFs or histograms")
    parser.add_argument('--no-fill', help="dry run -- do not fill histograms", action="store_true", 
            dest='noFill')
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    if tag not in ["Razor2015","Razor2016"]:
        sys.exit("Error: tag "+tag+" not supported!")

    #initialize
    sfHists = {}
    plotOpts = { 'comment':False }

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
        os.system('mkdir -p '+outdir)
        #get correct variable names
        sfVars = ("MR","Rsq")
        if region == "GJetsInv": sfVars = ("MR_NoPho","Rsq_NoPho")
        elif "Inv" in region: sfVars = ("MR_NoW", "Rsq_NoW")
        #QCD estimate for photon region
        if region == 'GJetsInv':
            dataDrivenQCD = True
        else:
            dataDrivenQCD = False
        #perform analysis
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, debugLevel=debugLevel, noFill=args.noFill, dataDrivenQCD=dataDrivenQCD) 
        #compute scale factors
        appendScaleFactors(process, hists, sfHists, lumiData=analysis.lumi, th2PolyXBins=xbins, 
                th2PolyCols=cols, debugLevel=debugLevel, var=sfVars, printdir=outdir)
        #export histograms
        if not args.noSave:
            macro.exportHists(hists, outFileName='controlHistograms'+region+'.root', 
                outDir=outdir, debugLevel=debugLevel)

    #write scale factors
    if not args.noSave:
        outfile = rt.TFile("data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), "RECREATE")
        for name in sfHists:
            print "Writing scale factor histogram",sfHists[name].GetName(),"to file"
            sfHists[name].Write(sfHists[name].GetName().replace("Poly",""))
        outfile.Close()
