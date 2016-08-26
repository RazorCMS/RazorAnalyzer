import sys,os,argparse,copy
import ROOT as rt

from macro import macro, razorWeights
from macro.razorAnalysis import Analysis
from macro.razorMacros import runFitAndToys, makeControlSampleHistsForAnalysis

FIT_DIR = "Fits"
TOYS_FILES = {
        "MultiJet":FIT_DIR+"/toys_Bayes_noStat_MultiJet.root",
        "LeptonMultiJet":FIT_DIR+"/toys_Bayes_noStat_LeptonMultiJet.root",
        "DiJet":FIT_DIR+"/toys_Bayes_noStat_DiJet.root",
        "LeptonJet":FIT_DIR+"/toys_Bayes_noStat_LeptonJet.root",
        }
FULL_TOYS_FILES = {
        }

commonShapeErrors = [
        ('singletopnorm',"SingleTop"),
        ('othernorm',"Other"),
        ('qcdnorm','QCD'),
        'btag', 'pileup', 'bmistag', 'facscale', 'renscale', 'facrenscale',
        ('btaginvcrosscheck',['ZInv']),
        ('btagcrosscheckrsq',['TTJets1L','TTJets2L','WJets']),
        ('btagcrosscheckmr',['TTJets1L','TTJets2L','WJets']),
        ('sfstatzinv',['ZInv']),
        ('sfsyszinv',['ZInv']),
        'jes','ees','mes',
        ('ttcrosscheck',['TTJets2L']),
        ('sfstatttjets',['TTJets1L','TTJets2L']),
        ('sfsysttjets',['TTJets1L','TTJets2L']),
        ('sfstatwjets',['WJets']),
        ('sfsyswjets',['WJets'])
        ]
lepShapeErrors = commonShapeErrors+['tightmuoneff','tighteleeff','muontrig','eletrig']
hadShapeErrors = commonShapeErrors+['vetolepptcrosscheck','vetotauptcrosscheck',
        'vetolepetacrosscheck','vetotauetacrosscheck','vetomuoneff','vetoeleeff']
shapes = { 'MultiJet':hadShapeErrors, 'LeptonMultiJet':lepShapeErrors, 
           'DiJet':hadShapeErrors, 'LeptonJet':lepShapeErrors, 
           }

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="display detailed output messages",
                                action="store_true")
    parser.add_argument("-d", "--debug", help="display excruciatingly detailed output messages",
                                action="store_true")
    parser.add_argument("--unblind", help="do not blind signal sensitive region", action='store_true')
    parser.add_argument('--no-mc', help="do not process MC, do data and fit only", 
            action='store_true', dest="noMC")
    parser.add_argument('--no-fit', help="do not load fit results, process data and MC only", 
            action='store_true', dest='noFit')
    parser.add_argument('--full', help="do full fit (default is sideband)", action='store_true')
    parser.add_argument('--no-data', help="do not process data, do fit and MC only", 
            action='store_true', dest='noData')
    parser.add_argument('--no-sys', help="no shape unncertainties or cross check systematics", 
            action="store_true", dest='noSys')
    parser.add_argument('--no-qcd', help="do not include QCD prediction", action="store_true", 
            dest='noQCD')
    parser.add_argument('--no-fill', help="dry run -- do not fill histograms", action="store_true", 
            dest='noFill')
    parser.add_argument('--no-sfs', help="ignore MC scale factors", action="store_true", 
            dest="noSFs")
    parser.add_argument('--box', help="choose a box")
    parser.add_argument('--btags', type=int, help="choose a number of btags")
    parser.add_argument('--b-inclusive', help='do not bin in btags', action='store_true',
            dest='bInclusive')
    parser.add_argument("--tag", dest="tag", required=True,
            help="Analysis tag, e.g. Razor2015")
    args = parser.parse_args()
    debugLevel = args.verbose + 2*args.debug
    tag = args.tag
    if tag not in ["Razor2015","Razor2016"]:
        sys.exit("Error: tag "+tag+" not supported!")

    #initialize
    plotOpts = {"SUS15004":True}

    doSideband=(not args.full)
    toysToUse = TOYS_FILES
    dirSuffix = ""
    if not doSideband:
        toysToUse = FULL_TOYS_FILES
        dirSuffix += 'Full'
        plotOpts['sideband'] = False
    else:
        plotOpts['sideband'] = True
    if not args.unblind:
        dirSuffix += 'Blinded'
    if args.noFit: 
        toysToUse = {}
        del plotOpts['sideband']
    boxesToUse = ["MultiJet", "LeptonMultiJet", "DiJet", "LeptonJet"]
    if args.box is not None:
        boxesToUse = [args.box]
    if args.btags is not None:
        btaglist = [args.btags]
    else:
        btaglist = [0,1,2,3]
    if args.noSFs:
        dirSuffix += 'NoSFs'
    if args.noSys:
        dirSuffix += 'NoSys'
    if args.bInclusive:
        dirSuffix += 'BInclusive'

    regionsOrder = []
    regions = {}
    for box in boxesToUse:
        for btags in btaglist:
            #deal with last b-tag bin
            if box in ["DiJet", "MuJet", "EleJet", "LeptonJet"]:
                if btags == 3: continue
                elif btags >= 2:
                    nbMax = -1 #no upper limit
                else:
                    nbMax = btags #exclusive
            else:
                if btags >= 3:
                    nbMax = -1 #no upper limit
                else:
                    nbMax = btags #exclusive
            if args.bInclusive:
                nbMin = 0
                nbMax = 10
            #define analysis region
            extBox = '%s%dB%s'%(box,btags,dirSuffix)
            regionsOrder.append(extBox)
            regions[extBox] = Analysis(box, tag=tag, nbMin=btags, nbMax=nbMax)

    ####LOAD ALL SCALE FACTOR HISTOGRAMS

    #scale factor file names
    sfdir = "data/ScaleFactors/RazorMADD2015/"
    sfFile = sfdir+'/RazorScaleFactors_%s.root'%(tag)
    sfFile_nJets = sfdir+'/RazorNJetsScaleFactors_%s.root'%(tag)
    vetolepFile = sfdir+'/RazorVetoLeptonClosureTests_%s.root'%(tag)
    ttFileDiJet = sfdir+'/TTBarDileptonSystematic_DiJet_%s.root'%(tag)
    ttFileMultiJet = sfdir+'/TTBarDileptonSystematic_MultiJet_%s.root'%(tag)
    dyFile = sfdir+'/RazorDYJetsDileptonInvCrossCheck_%s.root'%(tag)
    btagFile = sfdir+'/RazorBTagClosureTests_%s.root'%(tag)

    #get MR-Rsq scale factor histograms
    sfNames={
            "ZInv":"GJetsInv",
            "TTJets1L":"TTJets",
            "TTJets2L":"TTJets",
            }
    processNames = regions.itervalues().next().samples
    sfHists = macro.loadScaleFactorHists(sfFilename=sfFile, processNames=processNames, 
            scaleFactorNames=sfNames, debugLevel=debugLevel)
    #reopen the file and grab the ZNuNu up/down histograms
    #down scale factors are (gjets - (wjets-gjets))
    sfTFile = rt.TFile.Open(sfFile)
    sfHists['ZInvUp'] = sfTFile.Get('WJetsInvScaleFactors')
    sfHists['ZInvDown'] = sfTFile.Get('GJetsInvScaleFactors_Down') 
    #get njets scale factor histogram
    sfNJetsFile = rt.TFile.Open(sfFile_nJets)
    sfHists['NJets'] = sfNJetsFile.Get("NJetsCorrectionScaleFactors")
    sfHists['NJetsInv'] = sfNJetsFile.Get("GJetsScaleFactorVsNJets")
    #get veto lepton and tau scale factor histograms
    vlFile = rt.TFile.Open(vetolepFile)
    for ltype in ['VetoLepton','VetoTau']:
        for jtype in ['DiJet','MultiJet']:
            name = jtype+'For'+ltype
            sfHists[ltype+jtype+'PtUp'] = vlFile.Get(name+'ScaleFactors')
            sfHists[ltype+jtype+'PtDown'] = macro.invertHistogram(sfHists[ltype+jtype+'PtUp'])
            sfHists[ltype+jtype+'EtaUp'] = vlFile.Get(name+'PtCorrScaleFactors')
            sfHists[ltype+jtype+'EtaDown'] = macro.invertHistogram(sfHists[ltype+jtype+'EtaUp'])
    #get DYJets and TTBar Dilepton cross check scale factor histograms
    ttTFiles = { "DiJet":rt.TFile.Open(ttFileDiJet), "MultiJet":rt.TFile.Open(ttFileMultiJet) }
    for ttname, ttfile in ttTFiles.iteritems():
        histname = 'TTJetsDilepton'+ttname
        sfHists[histname+'Up'] = ttfile.Get('TTBarDileptonSystematic')
        #convert to correct SF histogram format
        for nb in range(sfHists[histname+'Up'].GetSize()+1):
            sfHists[histname+'Up'].SetBinContent( nb, 
                    sfHists[histname+'Up'].GetBinContent(nb)+1.0 )
        #get 'down' version of histogram
        sfHists[histname+'Down'] = macro.invertHistogram(sfHists[histname+'Up'])
    dyTFile = rt.TFile.Open(dyFile)
    for jtype in ['DiJet','MultiJet']:
        sfHists['DYJetsInv'+jtype+'Up'] = dyTFile.Get('DYJetsDileptonInv'+jtype+'ScaleFactors')
        sfHists['DYJetsInv'+jtype+'Down'] = macro.invertHistogram(sfHists['DYJetsInv'+jtype+'Up'])
    #get b-tag closure test results
    btagTFile = rt.TFile.Open(btagFile)
    for jtype in ['DiJet','MultiJet']:
        for b in range(4):
            if jtype == 'DiJet' and b > 2: continue
            bs = str(b)
            sfHists['MR'+jtype+bs+'BUp'] = btagTFile.Get('OneLepton'+jtype+'ClosureTest'+bs+'BMRScaleFactors')
            sfHists['Rsq'+jtype+bs+'BUp'] = btagTFile.Get('OneLepton'+jtype+'ClosureTest'+bs+'BRsqScaleFactors')
            sfHists['MR'+jtype+bs+'BDown'] = macro.invertHistogram(sfHists['MR'+jtype+bs+'BUp'])
            sfHists['Rsq'+jtype+bs+'BDown'] = macro.invertHistogram(sfHists['Rsq'+jtype+bs+'BUp'])
    #get ZInv b-tag cross check histogram
    sfHists['ZInvBUp'] = btagTFile.Get('ZNuNuBTagClosureSysUnc')
    #convert to correct SF histogram format
    for nb in range(sfHists['ZInvBUp'].GetSize()+1):
        sfHists['ZInvBUp'].SetBinContent( nb, sfHists['ZInvBUp'].GetBinContent(nb)+1.0 )
    #get 'down' version of histogram
    sfHists['ZInvBDown'] = macro.invertHistogram(sfHists['ZInvBUp'])

    #check that everything came out correctly
    for h,hist in sfHists.iteritems():
        if debugLevel > 0:
            print "Checking scale factor histogram:",h
        assert hist

    #estimate yields in signal region
    for region in regionsOrder:
        analysis = regions[region]
        boxName = region.replace(dirSuffix,'')[:-2]
        btags = int(region.replace(dirSuffix,'')[-2])
        #get correct NJets scale factors
        auxSFs = razorWeights.getNJetsSFs(analysis,jetName='nSelectedJets')

        print "\nBox:",region,"("+boxName,str(btags),"B-tag)"

        #make output directory
        outdir = "Plots/"+tag+"/"+region
        os.system('mkdir -p '+outdir)

        blindBins = [(x,y) for x in range(2,len(analysis.binning["MR"])+1) 
                for y in range(2,len(analysis.binning["Rsq"])+1)]

        #apply options
        if args.unblind: blindBins = None
        if args.noQCD and 'QCD' in analysis.samples:
            analysis.samples.remove('QCD')
        if args.noMC: analysis.samples = []
        if analysis.samples is None or len(analysis.samples) == 0:
            analysis.filenames = {"Data":analysis.filenames["Data"]}
        if args.noData: 
            del analysis.filenames['Data']
        shapesToUse = copy.copy(shapes[boxName])
        if args.noSys:
            shapesToUse = []
        if (boxName not in toysToUse or toysToUse[boxName] is None) and 'sideband' in plotOpts:
            del plotOpts['sideband']

        ##### get correct set of scale factor histograms
        sfHistsToUse = sfHists.copy()
        auxSFsToUse = auxSFs.copy()
        if 'MultiJet' in boxName:
            jtype = 'MultiJet'
        else:
            jtype = 'DiJet'
        #veto lepton/tau
        for ltype in ['VetoLepton','VetoTau']:
            for pteta in ['Pt','Eta']:
                for updown in ['Up','Down']:
                    sfHistsToUse[ltype+pteta+updown] = sfHistsToUse[ltype+jtype+pteta+updown]
        ##ttbar dilepton and dyjets dilepton
        for name in ['TTJetsDilepton','DYJetsInv']:
            for updown in ['Up','Down']:
                sfHistsToUse[name+updown] = sfHistsToUse[name+jtype+updown]
        #b-tag closure tests
        for b in range(4):
            bs = str(b)
            if jtype == 'DiJet' and b > 2: continue
            for updown in ['BUp','BDown']:
                for mrrsq in ['MR','Rsq']:
                    sfHistsToUse[mrrsq+bs+updown] = sfHistsToUse[mrrsq+jtype+bs+updown]

        #option to disable scale factors
        if args.noSFs:
            print "Ignoring all scale factor histograms and uncertainties from scale factor cross checks."
            sfHistsToUse = {}
            auxSFsToUse = {}
            toRemove = ['btaginvcrosscheck','btagcrosscheckrsq','btagcrosscheckmr','sfsyszinv','ttcrosscheck','zllcrosscheck','sfsysttjets','sfsyswjets','vetolepptcrosscheck','vetotauptcrosscheck','vetolepetacrosscheck','vetotauetacrosscheck']
            #remove scale factor cross check uncertainties
            shapesToUse = [s for s in shapesToUse if s not in toRemove]
            #this removes scale factor uncertainties that are listed as tuples
            shapesToUse = [s for s in shapesToUse if not (hasattr(s, '__getitem__') and s[0] in toRemove)] 

        #run analysis
        hists = makeControlSampleHistsForAnalysis( analysis,
                sfHists=sfHistsToUse, treeName="RazorInclusive", 
                shapeErrors=shapesToUse, fitToyFiles=toysToUse, boxName=boxName, blindBins=blindBins,
                btags=btags, debugLevel=debugLevel, auxSFs=auxSFsToUse, dataDrivenQCD=True, printdir=outdir, 
                plotOpts=plotOpts, noFill=args.noFill, exportShapeErrs=True, propagateScaleFactorErrs=False)
        #export histograms
        macro.exportHists(hists, outFileName='razorHistograms'+region+'.root', outDir=outdir, 
                debugLevel=debugLevel)
