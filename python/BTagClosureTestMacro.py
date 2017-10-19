import sys, os
import ROOT as rt

import macro.razorAnalysis as razor
from macro import macro, razorWeights
from macro.razorMacros import makeControlSampleHistsForAnalysis, appendScaleFactors

def getOutputFilename(tag='Razor2016_MoriondRereco'):
    return "data/ScaleFactors/RazorMADD2015/RazorBTagScaleFactors_%s.root"%(tag)

def loadScaleFactors(sfHists={}, tag='Razor2016_MoriondRereco', gjets=False):
    f = rt.TFile.Open(getOutputFilename(tag))
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
    for jets in ['DiJet', 'MultiJet']:
        for btags in range(4):
            if (jets == 'DiJet' or gjets) and btags > 2:
                continue
            name = jets+str(btags)+'BScaleFactors'
            h = f.Get(sfKey+name)
            assert(h)
            h.SetDirectory(0)
            sfHists[sfKey+jets+str(btags)+'B'] = h
    return sfHists

def isMultiJet(analysis):
    return analysis.njetsMin >= 4 or 'MultiJet' in analysis.region

def getSFHistName(analysis, gjets=False):
    jetName = 'DiJet'
    if isMultiJet(analysis):
        jetName = 'MultiJet'
    nbtags = max(0, analysis.nbMin)
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
        if nbtags > 2:
            nbtags = 2
    return sfKey+jetName+str(nbtags)+'B'

def adjustForRegion(analysis, sfHists, auxSFs, gjets=False):
    """
    Adjusts scale factor dictionaries so that the correct
    corrections are applied for the current region
    """
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
    if sfKey in sfHists:
        # This avoids using last region's scale factors
        del sfHists[sfKey]
    sfHistName = getSFHistName(analysis, gjets=gjets)
    if sfHistName in sfHists:
        # This sets up the correct reweighting histogram
        sfHists[sfKey] = sfHists[sfHistName]
        print "Using histogram {} for {}".format(sfHistName, sfKey)
    else:
        # This avoids trying to apply SFs we don't have
        for proc, sfs in auxSFs.iteritems():
            if sfKey in sfs:
                print "Omitting {} scale factors for {}".format(sfKey, proc)
                del sfs[sfKey]

def adjustForRegionBInclusive(analysis, sfHists, auxSFs, gjets=False):
    """
    Same as adjustForRegion, but sets up all b-tag 
    histograms simulataneously.  Use to process a sample
    that contains events of different b-tag multiplicities.
    """
    sfKey = 'MR'
    if gjets:
        sfKey = 'MRInv'
    nbMax = 2
    jets = 'DiJet'
    if isMultiJet(analysis):
        if not gjets:
            nbMax = 3
        jets = 'MultiJet'
    for nb in range(nbMax + 1):
        thisSFKey = sfKey+str(nb)+'B'
        if thisSFKey in sfHists:
            # This avoids using last region's scale factors
            del sfHists[thisSFKey]
        sfHistName = sfKey+jets+str(nb)+'B'
        if sfHistName in sfHists:
            # This sets up the correct reweighting histogram
            sfHists[thisSFKey] = sfHists[sfHistName]
            print "Using histogram {} for {}".format(sfHistName, thisSFKey)
        else:
            # This avoids trying to apply SFs we don't have
            for proc, sfs in auxSFs.iteritems():
                if thisSFKey in sfs:
                    print "Omitting {} scale factors for {}".format(thisSFKey, proc)
                    del sfs[thisSFKey]


if __name__ == "__main__":
    rt.gROOT.SetBatch()

    parser = razor.make_parser()
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    boostCuts = not args.noBoostCuts

    plotOpts = { "comment":False, 'SUS15004CR':True } 
    regions = {}
    regionsOrder = []

    #define all tests
    jetsOrder = ["DiJet","MultiJet"]
    jetsLimit = [(2,3),(4,-1)]
    for name,jets in zip(jetsOrder, jetsLimit):
        for nb in range(4):
            if nb >= 3 and name == 'DiJet': 
                continue
            regionName = 'OneLepton'+name+str(nb)+'B'
            regions[regionName] = razor.Analysis("SingleLepton",tag=tag,
                    njetsMin=jets[0], njetsMax=jets[1], nbMin=nb, nbMax=nb, 
                    boostCuts=boostCuts)
            regionsOrder.append(regionName)
            regions[regionName+'MRCorr'] = razor.Analysis("SingleLepton",tag=tag,
                    njetsMin=jets[0], njetsMax=jets[1], nbMin=nb, nbMax=nb, 
                    boostCuts=boostCuts)
            regionsOrder.append(regionName+'MRCorr')

    sfHists = macro.loadScaleFactorHists(
            sfFilename="data/ScaleFactors/RazorMADD2015/RazorScaleFactors_%s.root"%(tag), 
            processNames=regions.itervalues().next().samples, debugLevel=debugLevel)
    sfVars = ("MR","Rsq")
    sfNJetsFile = rt.TFile.Open(
            "data/ScaleFactors/RazorMADD2015/RazorNJetsScaleFactors_%s.root"%(tag))
    sfHists['NJetsTTJets'] = sfNJetsFile.Get("TTJetsScaleFactors")
    sfHists['NJetsWJets'] = sfNJetsFile.Get("WJetsScaleFactors")

    for region in regionsOrder:
        print "\nRegion:",region,"\n"
        outdir = 'Plots/'+tag+'/'+region
        os.system('mkdir -p '+outdir)
        analysis = regions[region]
        auxSFs = razorWeights.getNJetsSFs(analysis)
        auxSFs = razorWeights.addBTagSFs(analysis, auxSFs)
        adjustForRegion(analysis, sfHists, auxSFs)
        hists = makeControlSampleHistsForAnalysis( analysis, plotOpts=plotOpts, sfHists=sfHists,
                sfVars=sfVars, printdir=outdir, auxSFs=auxSFs, btags=analysis.nbMin,
                debugLevel=debugLevel, noFill=args.noFill )

        sfHistName = getSFHistName(analysis)
        if 'MRCorr' in region:
            sfHistName = sfHistName.replace('MR', 'Rsq')
            appendScaleFactors( sfHistName, hists, sfHists, lumiData=analysis.lumi, var="Rsq",
                    debugLevel=debugLevel, signifThreshold=1.0, printdir=outdir )
        else:
            appendScaleFactors( sfHistName, hists, sfHists, lumiData=analysis.lumi, var="MR",
                    debugLevel=debugLevel, printdir=outdir )
        if not args.noSave:
            macro.exportHists( hists, outFileName='controlHistograms'+region+'.root',
                    outDir=outdir, debugLevel=debugLevel )
            outfile = rt.TFile(getOutputFilename(tag), "UPDATE")
            histToWrite = sfHists[sfHistName]
            print "Writing scale factor histogram",histToWrite.GetName(),"to file"
            outfile.cd()
            histToWrite.Write( histToWrite.GetName() )
            outfile.Close()
