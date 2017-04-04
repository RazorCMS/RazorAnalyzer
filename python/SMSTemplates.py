import sys
import argparse
import copy
import ROOT as rt

#local imports
from DustinTuples2DataCard import uncorrelate, getTheoryCrossSectionAndError
from macro.razorAnalysis import Analysis
from macro.razorMacros import makeControlSampleHistsForAnalysis, getMaxBtags, unrollAndStitch

signalShapeUncerts = ['tightmuoneff','tighteleeff','vetomuoneff',
        'vetoeleeff','jes','muontrig','eletrig','btag',
        'tightmuonfastsim','tightelefastsim','vetomuonfastsim',
        'vetoelefastsim','btagfastsim','facscale','renscale',
        'facrenscale','ees','mes','pileup','isr']

def getModelInfoFromFilename(f, xBR=-1, yBR=-1):
    """f: filename
       xBR, yBR: branching fractions
       returns: model_name, mass1, mass2"""
    modelString = '_'.join(f.split('/')[-1].split(
        '.root')[0].split('_')[:-2])
    if xBR>-1 and yBR>-1:
        modelString = modelString.replace('T1ttbb',
                ('T1x%.2fy%.2f'%(xBR,yBR)).replace('.','p'))
    modelName = modelString.split('-')[-1]
    massPoint = '_'.join(f.split('/')[-1].split(
        '.root')[0].split('_')[1:])
    m1 = massPoint.split('_')[-2]
    m2 = massPoint.split('_')[-1]
    return modelName, m1, m2

def makeSMSTemplates(box, inFile, uncertainties=[], debugLevel=0,
        xBR=-1, yBR=-1, tag="Razor2016_MoriondRereco"):
    """Returns a dictionary of histograms representing predicted
        yields for the indicated signal sample"""
    print '\nInput file is %s' % inFile
    uncerts = copy.copy(uncertainties)
    doMCStat = False
    if 'mcstat%s'%box.lower() in uncerts:
        doMCStat = True
        uncerts.remove('mcstat%s'%box.lower())

    # get mass point information
    model, m1, m2 = getModelInfoFromFilename(inFile, xBR, yBR)
    isGluinos = ("T1" in model or "T5" in model)
    if isGluinos:
        thyXsec, _ = getTheoryCrossSectionAndError(
                mGluino=m1)
    else:
        thyXsec, _ = getTheoryCrossSectionAndError(
                mStop=m1)
    print "\n--- %s %s %s ---\n"%(model, m1, m2)
    print "Theory cross section: %.3f pb"%thyXsec

    minBtags = 0
    maxBtags = getMaxBtags(box)
    # special case: 1L control regions
    if box == "TTJetsSingleLeptonForSignal":
        minBtags = maxBtags = 1
    elif box == "WJetsSingleLeptonForSignal":
        minBtags = maxBtags = 0
    hists = []
    unrollBins = []
    for nb in range(minBtags, maxBtags+1):
        # get analysis info
        nbMax = nb
        if nb == maxBtags:
            nbMax = -1
        # special case: 1L control regions
        if box == "WJetsSingleLeptonForSignal":
            nbMax = 0
        analysis = Analysis(box, tag, nbMin=nb, nbMax=nbMax)
        unrollBins.append(analysis.unrollBins)

        # modify for signal sample
        analysis.filenames = { 'Signal':inFile }
        analysis.samples = ['Signal']
        analysis.samplesReduced = analysis.samples
            
        # scale according to cross section
        f = rt.TFile.Open(inFile)
        nEvents = f.Get('NEvents').Integral()
        globalScaleFactor = thyXsec/nEvents 
        print "Number of events: %d"%nEvents
        print "Integrated luminosity: %d /pb"%analysis.lumi
        print "Overall scale factor: %.3f"%(
                analysis.lumi * globalScaleFactor)

        # fill the histograms
        hists.append(makeControlSampleHistsForAnalysis(analysis, 
                treeName="RazorInclusive", shapeErrors=uncerts,
                boxName=box, btags=nb, makePlots=False, 
                exportShapeErrs=True, propagateScaleFactorErrs=False,
                lumiMC=1./globalScaleFactor, debugLevel=debugLevel))
            
    signalHists = unrollAndStitch(hists, box, samples=analysis.samples,
            debugLevel=debugLevel, unrollBins=unrollBins)

    #perform uncorrelation procedure (for MC stat uncertainty)
    if doMCStat:
        uncorrelate(signalHists, 'mcstat%s'%box.lower())

    return signalHists

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inFile', help='Input file path')
    parser.add_argument('--tag', default='Razor2016_MoriondRereco')
    parser.add_argument('-b','--box', default="MultiJet",
                  help="box name")
    parser.add_argument('--xBR', default=-1, type=float,
                  help="x = BR(~g -> b b ~chi0)")
    parser.add_argument('--yBR', default=-1,type=float,
                  help="y = BR(~g -> t t ~chi0)")
    parser.add_argument('--no-signal-sys', dest="noSignalSys", 
                  action='store_true', 
                  help="no signal systematic templates")
    args = parser.parse_args()
    
    boxList = box.split('_')

    if args.noSignalSys:
        shapes = []
    else:
        shapes = signalShapeUncerts

    for curBox in boxList:
        print "Testing SMS template maker for box", args.box
        hists = makeSMSTemplates(args.box, args.inFile, shapes,
                args.xBR, args.yBR, args.tag)
