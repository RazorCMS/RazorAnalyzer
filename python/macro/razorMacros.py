## Machinery for inclusive razor analysis

import sys,os
from array import array
import ROOT as rt
from collections import namedtuple

#local imports
from RunCombine import exec_me
from framework import Config
import macro
import plotting
from razorAnalysis import *
from razorWeights import *

###########################################
### RAZOR FIT
###########################################

def runFitAndToys(fitDir, boxName, lumi, dataName, dataDir='./', config='config/run2.config', sideband=False, numToys=4000, noStat=False):
    #make folder
    if not os.path.isdir(fitDir):
        exec_me('mkdir -p '+fitDir, False)
    #make RooDataSet
    exec_me('python python/DustinTuple2RooDataSet.py -b '+boxName+' -c '+config+' -l '+str(lumi)+' --data -d '+fitDir+' '+dataDir+'/'+dataName+'.root', False)
    #do fit
    if not sideband:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' --data '+fitDir+'/'+dataName+'_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    else:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' --data --fit-region LowMR,LowRsq '+fitDir+'/'+dataName+'_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    #run toys
    exec_me('python python/RunToys.py -b '+boxName+' -c '+config+' -i '+fitDir+'/BinnedFitResults_'+boxName+'.root -d '+fitDir+' -t '+str(numToys)+((noStat)*" --no-stat"), False)

def runFitAndToysMC(fitDir, boxName, lumi, fileNames, mcDir='./', config='config/run2.config', sideband=False, numToys=4000, noStat=False):
    #make folder
    if not os.path.isdir(fitDir):
        exec_me('mkdir -p '+fitDir, False)
    #make RooDataSet
    exec_me('python python/DustinTuple2RooDataSet.py -w -b '+boxName+' -c '+config+' -l '+str(lumi)+' -d '+fitDir+' '+' '.join(fileNames), False)
    #do fit
    if not sideband:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' '+fitDir+'/RazorInclusive_SMCocktail_weighted_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    else:
        exec_me('python python/BinnedFit.py -c '+config+' -d '+fitDir+' -l '+str(lumi)+' -b '+boxName+' --fit-region LowMR,LowRsq '+fitDir+'/RazorInclusive_SMCocktail_weighted_lumi-'+('%1.3f' % (lumi*1.0/1000))+'_0-3btag_'+boxName+'.root', False)
    #run toys
    exec_me('python python/RunToys.py -b '+boxName+' -c '+config+' -i '+fitDir+'/BinnedFitResults_'+boxName+'.root -d '+fitDir+' -t '+str(numToys)+((noStat)*" --no-stat"), False)

def get2DNSigmaHistogram(data, bins, fitToyFiles, boxName, btags=-1, debugLevel=0):
    print "Making Nsigma histogram using fit information"
    #set up histogram
    nsigma = data.Clone(data.GetName()+"Nsigma")
    nsigma.Reset()

    #load fit information, including toys
    toyFile = rt.TFile(fitToyFiles[boxName])
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFiles[boxName],"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    #make dummy options tuple
    Opt = namedtuple("Opt", ["printErrors","noStat"])
    options = Opt(False,False)

    z = [0.,1.,2.,3.,4.] #b-tag binning
    if btags < 0: #be inclusive in b-tags
        zmin = 0
        zmax = len(z)-1
    else:
        zmin = btags+1
        zmax = btags+1

    from PlotFit import getNsigma2D
    nsigma = getNsigma2D(nsigma, data, toyTree, options, 'yx', 0, len(bins['MR'])-1, 0, len(bins['Rsq'])-1, zmin, zmax, bins['MR'], bins['Rsq'], z)
    return nsigma

def getFitCorrelationMatrix(config, box, fitToyFile, debugLevel=0):
    print "Creating fit correlation matrix"

    cfg = Config.Config(config)
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning

    #load fit information, including toys
    toyFile = rt.TFile(fitToyFile)
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFile,"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    from PlotFit import getCorrelationMatrix
    corrMatrix = getCorrelationMatrix(toyTree,'yx',0,len(x)-1,0,len(y)-1,0,len(z)-1,x,y,z)
    return corrMatrix

def import2DRazorFitHistograms(hists, bins, fitToyFile, c, dataName="Data", btags=-1, debugLevel=0, noStat=False):
    print "Loading fit histograms"
    if noStat: 
        print "Using only systematic errors on fit points"
        if 'noStat' not in fitToyFile:
            fitToyFile = fitToyFile.replace('Bayes','Bayes_noStat')
    #sanity check
    if "Fit" in hists:
        print "Error in import2DFitHistograms: fit histogram already exists!"
        return
    #make histogram for fit result
    hists["Fit"] = {}
    for v in ["MR","Rsq",("MR","Rsq")]:
        hists["Fit"][v] = next(hists.itervalues())[v].Clone(next(hists.itervalues())[v].GetName()+"Fit")
    hists["Fit"][v].Reset()

    #load fit information, including toys
    toyFile = rt.TFile(fitToyFile)
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFile,"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    #get uncertainties on fit prediction
    #using code imported from Javier's plotting script
    from PlotFit import getBinSumDicts, getBestFitRms, getErrors1D
    z = [0.,1.,2.,3.,4.] #b-tag binning
    if btags < 0: #be inclusive in b-tags
        zmin = 0
        zmax = len(z)-1
    else:
        zmin = btags+1
        zmax = btags+1

    #make dummy options tuple
    Opt = namedtuple("Opt", ["printErrors","noStat"])
    options = Opt(False, noStat)

    #store best fit and uncertainty in each bin
    #1D
    binSumDict = getBinSumDicts('x', 0,len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
    for i,sumName in binSumDict.iteritems():
        if dataName in hists: nObs = hists[dataName]["MR"].GetBinContent(i)
        else: nObs = 0
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options,"")
        hists["Fit"]["MR"].SetBinContent(i,bestFit)
        hists["Fit"]["MR"].SetBinError(i,rms)
    binSumDict = getBinSumDicts('y', 0, len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
    for i,sumName in binSumDict.iteritems():
        if dataName in hists: nObs = hists[dataName]["Rsq"].GetBinContent(i)
        else: nObs = 0
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options,"")
        hists["Fit"]["Rsq"].SetBinContent(i,bestFit)
        hists["Fit"]["Rsq"].SetBinError(i,rms)
    #2D
    binSumDict = getBinSumDicts('yx', 0, len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
    for (i,j),sumName in binSumDict.iteritems():
        if dataName in hists: nObs = hists[dataName][("MR","Rsq")].GetBinContent(i,j)
        else: nObs = 0
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options=options,plotName="")
        hists["Fit"][("MR","Rsq")].SetBinContent(i,j,bestFit)
        hists["Fit"][("MR","Rsq")].SetBinError(i,j,rms)

def get3DRazorFitHistogram(configFile, fitToyFile, boxName, debugLevel=0):
    print "Getting 3D fit histogram from",fitToyFile
    #load binning
    cfg = Config.Config(configFile)
    binsX = array('d', cfg.getBinning(boxName)[0])
    binsY = array('d', cfg.getBinning(boxName)[1])
    binsZ = array('d', cfg.getBinning(boxName)[2])
    nBins = (len(binsX)-1)*(len(binsY)-1)*(len(binsZ)-1)

    #load functions from Javier
    from PlotFit import getBinSumDicts, getBestFitRms

    #make fit histogram 
    fitHist3D = rt.TH3F("fit"+boxName, "fit"+boxName, len(binsX)-1, binsX, len(binsY)-1, binsY, len(binsZ)-1, binsZ)

    #load fit information, including toys
    toyFile = rt.TFile(fitToyFile)
    assert toyFile
    if debugLevel > 0: print "Opened file",fitToyFile,"to get Bayesian toy results"
    toyTree = toyFile.Get("myTree")
    assert toyTree
    if debugLevel > 0: print "Got tree myTree"

    #make dummy options tuple
    Opt = namedtuple("Opt", ["printErrors","noStat"])
    options = Opt(False, False)

    #store best fit and uncertainty in each bin
    binSumDict = getBinSumDicts('zyx', 0,len(binsX)-1,0,len(binsY)-1,0,len(binsZ)-1,binsX,binsY,binsZ)
    c = rt.TCanvas("placeholder","placeholder",800,600)
    for (i,j,k),sumName in binSumDict.iteritems():
        bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,0,c,options,"")
        fitHist3D.SetBinContent(i,j,k,bestFit)
        fitHist3D.SetBinError(i,j,k,rms)

    return fitHist3D

def makeRazor3DTable(hist, boxName, signalHist=None, signalName="T1bbbb"):
    """Print latex table with prediction and uncertainty in each bin"""
    xbinLowEdges = []
    xbinUpEdges = []
    ybinLowEdges = []
    ybinUpEdges = []
    zbinLowEdges = []
    zbinUpEdges = []
    predictions = []
    uncerts = []
    signal = []
    #for each bin, get values for all table columns
    for bx in range(1, hist.GetNbinsX()+1):
        for by in range(1, hist.GetNbinsY()+1):
            for bz in range(1, hist.GetNbinsZ()+1):
                xbinLowEdges.append('%.0f' % (hist.GetXaxis().GetBinLowEdge(bx)))
                xbinUpEdges.append('%.0f' % (hist.GetXaxis().GetBinUpEdge(bx)))
                ybinLowEdges.append(str(hist.GetYaxis().GetBinLowEdge(by)))
                ybinUpEdges.append(str(hist.GetYaxis().GetBinUpEdge(by)))
                zbinLowEdges.append('%.0f' % (hist.GetZaxis().GetBinLowEdge(bz)))
                zbinUpEdges.append('%.0f' % (hist.GetZaxis().GetBinUpEdge(bz)))
                pred = hist.GetBinContent(bx,by,bz)
                predictions.append('%.2f' % (pred))
                unc = hist.GetBinError(bx,by,bz)
                if pred != 0:
                    uncerts.append('%.1f\\%%' % (100*(unc/pred)))
                else:
                    uncerts.append('%.2f' % (unc))
                if signalHist is not None:
                    signal.append('%.2f' % (signalHist.GetBinContent(bx,by,bz)))
    xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
    yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
    zRanges = copy.copy(zbinLowEdges)
    headers=["$M_R$", "$R^2$", "B-tags", "Prediction", "Uncertainty"]
    cols = [xRanges, yRanges, zRanges, predictions, uncerts]
    if signalHist is not None:
        headers.append(signalName)
        cols.append(signal)
    plotting.table_basic(headers, cols, caption="Fit prediction for the "+boxName+" box", printstr="razorFitTable"+boxName)

def makeRazor2DTable(pred, obs, nsigma, boxName, btags=-1):
    """Print latex table with prediction and uncertainty in each bin"""
    xbinLowEdges = []
    xbinUpEdges = []
    ybinLowEdges = []
    ybinUpEdges = []
    zbinLowEdges = []
    preds = []
    uncerts = []
    obses = []
    nsigmas = []
    #for each bin, get values for all table columns
    for bx in range(1, pred.GetNbinsX()+1):
        for by in range(1, pred.GetNbinsY()+1):
            xbinLowEdges.append('%.0f' % (pred.GetXaxis().GetBinLowEdge(bx)))
            xbinUpEdges.append('%.0f' % (pred.GetXaxis().GetBinUpEdge(bx)))
            ybinLowEdges.append(str(pred.GetYaxis().GetBinLowEdge(by)))
            ybinUpEdges.append(str(pred.GetYaxis().GetBinUpEdge(by)))
            zbinLowEdges.append('%.0f' % (btags))
            prediction = pred.GetBinContent(bx,by)
            uncert = pred.GetBinError(bx,by)
            preds.append('%.2f' % (prediction))
            uncerts.append('%.2f' % (uncert))
            observed = obs.GetBinContent(bx,by)
            obses.append('%.2f' % (observed))
            nsig = nsigma.GetBinContent(bx,by)
            nsigmas.append('%.2f' % (nsig))
    xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
    yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
    zRanges = copy.copy(zbinLowEdges)
    caption = "Comparison of observed and expected event yields for the "+boxName+" box"
    if btags >= 0:
        headers=["$M_R$", "$R^2$", "B-tags", "Prediction", "Uncertainty", "Observed", "Number of sigmas"]
        cols = [xRanges, yRanges, zRanges, preds, uncerts, obses, nsigmas]
        caption += " ("+str(btags)+" b-tags)"
    else:
        headers=["$M_R$", "$R^2$", "Prediction", "Uncertainty", "Observed", "Number of sigmas"]
        cols = [xRanges, yRanges, preds, uncerts, obses, nsigmas]
    plotting.table_basic(headers, cols, caption=caption, printstr="razor2DFitTable"+boxName+str(btags)+"btag")

###########################################
### BASIC HISTOGRAM FILLING/PLOTTING MACRO
###########################################

def makeControlSampleHists(regionName="TTJetsSingleLepton", filenames={}, samples=[], cutsMC="", cutsData="", bins={}, plotOpts={}, lumiMC=1, lumiData=3000, weightHists={}, sfHists={}, treeName="ControlSampleEvent",dataName="Data", weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], shapeErrors=[], miscErrors=[], fitToyFiles=None, boxName=None, btags=-1, blindBins=None, makePlots=True, debugLevel=0, printdir=".", plotDensity=True, sfVars = ("MR","Rsq")):
    """Basic function for filling histograms and making plots.

    regionName: name of the box/bin/control region (used for plot labels)
    filenames: dictionary of process:filename pairs for loading ntuples
    samples: list of samples to process, in the order that they should appear in stacked plots, legends, etc
    cutsMC, cutsData: strings, to be used with TTreeFormula to make selection cuts
    bins: dictionary formatted like { "variable1":[bin0,bin1,bin2,...], "variable2":[bin0,bin1,bin2,...]}
    plotOpts: optional -- dictionary of misc plotting options (see below for supported options)
    lumiMC, lumiData: in /pb
    weightHists: dictionary of weight histograms, like that produced by razorWeights.loadWeightHists
    sfHists: dictionary of scale factor histograms, like that produced by razorWeights.loadScaleFactorHists
    treeName: name of the tree containing input events
    weightOpts: list of strings with directives for applying weights to the MC
    shapeErrors: list of MC shape uncertainties
    miscErrors: optional -- list of misc uncertainty options (see below for supported options)
    fitToyFiles: optional -- dictionary of boxName:toyFile pairs for loading razor fit results
    boxName: optional -- name of razor box
    btags: optional -- used only to specify which fit to load"""

    titles = {
        "MR": "M_{R} (GeV)", 
        "Rsq": "R^{2}",
        "mll": "m_{ll} (GeV)",
        "NBJetsMedium" : "Number of B-tagged Jets",
        "NJets80" : "Number of Jets with p_{T} > 80 GeV",
        "NJets40" : "Number of Jets",        
        "lep1.Pt()": "lepton p_{T} (GeV)",
        "lep2.Pt()": "lepton p_{T} (GeV)",
        "lep1.Eta()": "lepton #eta",
        "lep2.Eta()": "lepton #eta",
        }

    ##Get plotting options (for customizing plot behavior)
    special = ""
    #set log scale
    if "logx" in plotOpts: logx = plotOpts["logx"]
    else: logx = True
    if "ymin" in plotOpts: ymin = plotOpts["ymin"]
    else: ymin = 0.1
    #allow disabling comment string (normally written at the top of each plot)
    if "comment" in plotOpts: comment = plotOpts["comment"]
    else: comment = True
    #use sideband fit result 
    if "sideband" in plotOpts:
        if plotOpts['sideband']:
            special += 'sideband'
        else:
            special += 'full'

    #set up files and trees
    inputs = filenames
    files = {name:rt.TFile.Open(inputs[name]) for name in inputs} #get input files
    for name in inputs: 
        assert files[name] #check input files
        if debugLevel > 0: print "Opened file",inputs[name]
    trees = macro.makeTreeDict(files, treeName, debugLevel)

    #add noise filters to cuts if bits are present in tree
    if dataName in trees: cutsData = appendNoiseFilters(cutsData, trees[dataName]) 

    #define histograms to fill
    hists,shapeHists = macro.setupHistograms(regionName, inputs, samples, bins, titles, shapeErrors, dataName)
    listOfVars = hists.itervalues().next().keys() #list of the variable names
    
    #fill histograms by looping over all trees
    if dataName in trees:
        print("\nData:")
        macro.loopTree(trees[dataName], weightF=weight_data, cuts=cutsData, hists=hists[dataName], weightHists=weightHists, weightOpts=[], debugLevel=debugLevel) 

    print("\nMC:")
    macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:hists[name] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, sfVars=sfVars, sysVars=sfVars, debugLevel=debugLevel) 

    #get up/down histogram variations
    for shape in shapeErrors:
        print "\n"+shape,"Up:"
        macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:shapeHists[name][shape+"Up"] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, errorOpt=shape+"Up", boxName=boxName, sfVars=sfVars, sysVars=None, debugLevel=debugLevel)
        print "\n"+shape,"Down:"
        macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:shapeHists[name][shape+"Down"] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, errorOpt=shape+"Down", boxName=boxName, sfVars=sfVars, sysVars=None, debugLevel=debugLevel)

    #propagate up/down systematics to central histograms
    macro.propagateShapeSystematics(hists, samples, bins, shapeHists, shapeErrors, miscErrors, boxName, debugLevel=debugLevel)

    c = rt.TCanvas(regionName+"c", regionName+"c", 800, 600)

    #optionally compare with fit
    nsigmaFitData = None
    nsigmaFitMC = None
    if fitToyFiles and "MR" in bins and "Rsq" in bins:
        import2DRazorFitHistograms(hists, bins, fitToyFiles[boxName], c, dataName, btags, debugLevel, noStat=True)
        if dataName in hists: 
            nsigmaFitData = get2DNSigmaHistogram(hists[dataName][("MR","Rsq")], bins, fitToyFiles, boxName, btags, debugLevel)
            print "Making nsigma histogram using data and fit prediction"
            makeRazor2DTable(pred=hists["Fit"][("MR","Rsq")], obs=hists[dataName][("MR","Rsq")],
                    nsigma=nsigmaFitData, boxName=boxName, btags=btags)

        if len(samples) > 0: #compare fit with MC
            #make total MC histogram
            mcTotal = hists[samples[0]][("MR","Rsq")].Clone("mcTotal")
            mcTotal.Reset()
            for name in samples:
                mcTotal.Add(hists[name][("MR","Rsq")])
            nsigmaFitMC = get2DNSigmaHistogram(mcTotal, bins, fitToyFiles, boxName, btags, debugLevel)
            print "Making nsigma histogram using MC and fit prediction"

    #print histograms
    rt.SetOwnership(c, False)
    if makePlots: macro.basicPrint(hists, mcNames=samples, varList=listOfVars, c=c, printName=regionName, logx=logx, dataName=dataName, ymin=ymin, comment=comment, lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", boxName=boxName, btags=btags, blindBins=blindBins, nsigmaFitData=nsigmaFitData, nsigmaFitMC=nsigmaFitMC, printdir=printdir, doDensity=plotDensity, special=special, vartitles=titles)

    #close files and return
    for f in files: files[f].Close()
    return hists

#######################################
### MAKE SCALE FACTORS FROM HISTOGRAMS
#######################################

def appendScaleFactors(process="TTJets", hists={}, sfHists={}, var=("MR","Rsq"), dataName="Data", normErrFraction=0.2, printTable=True, lumiData=0, signifThreshold=0., debugLevel=0, printdir="."):
    """Subtract backgrounds and make the data/MC histogram for the given process.
    Also makes up/down histograms corresponding to uncertainty on the background normalization (controlled by the normErrFraction argument).
    
    process: MC physics process for which scale factors should be computed
    hists: dictionary of data and MC histograms like that produced by the macro.loopTrees function
    sfHists: dictionary of existing scale factor histograms. the new scale factor histograms will be inserted into this dictionary.
    var: variable or tuple of variables in which scale factors should be computed (usually ("MR","Rsq") is used)
    signifThreshold: scale factors that are within N sigma of 1.0, where N=signifThreshold, are set to 1.
    """

    ##Sanity checks
    if debugLevel > 0: 
        print "Scale factor histograms so far:"
        print sfHists
    #warn if this scale factor histogram already exists
    if process in sfHists:
        print "Warning in appendScaleFactors: ",process," scale factor histogram already exists!  Will overwrite..."
    #warn if the needed histograms are not found
    if process not in hists:
        print "Error in appendScaleFactors: target MC histogram (",process,") was not found!"
        return
    if var not in hists[process]:
        print "Error in appendScaleFactors: could not find ",var," in hists[",process,"]!"
        return
    if dataName not in hists:
        print "Error in appendScaleFactors: target data histogram (",dataName,") was not found!"
        return
    if var not in hists[dataName]:
        print "Error in appendScaleFactors: could not find ",var," in hists[",dataName,"]!"


    #make the scale factor histogram (clone data hist; later subtract backgrounds and divide by MC hist)
    print "Making scale factor histogram for",process
    sfHists[process] = hists[dataName][var].Clone(process+"ScaleFactors")
    sfHists[process].SetDirectory(0)
    #make up/down systematic error histograms (this is the systematic due to MC background normalization)
    sfHists[process+"NormUp"] = hists[dataName][var].Clone(process+"ScaleFactorsUp")
    sfHists[process+"NormUp"].SetDirectory(0)
    sfHists[process+"NormDown"] = hists[dataName][var].Clone(process+"ScaleFactorsDown")
    sfHists[process+"NormDown"].SetDirectory(0)

    #subtract backgrounds in data
    bgProcesses = [mcProcess for mcProcess in hists if mcProcess != process and mcProcess != dataName and mcProcess != "Fit"] #get list of non-data, non-fit, non-signal samples
    for mcProcess in bgProcesses:

        #make sure relevant background histogram exists
        if var not in hists[mcProcess]:
            print "Error in appendScaleFactors: could not find",var," in hists[",mcProcess,"]!  Returning from appendScaleFactors..."
            return
        #subtract it
        if debugLevel > 0: print "Subtracting",mcProcess,"from",dataName,"distribution"
        sfHists[process].Add(hists[mcProcess][var], -1) 
        if mcProcess not in sfHists: #if we have not computed scale factors for this process, apply a flat normalization uncertainty to its yield
            if debugLevel > 0: print "Process",mcProcess,"has no associated scale factors.  Its normalization will be given a ",int(normErrFraction*100),"% uncertainty"
            sfHists[process+"NormUp"].Add(hists[mcProcess][var], -(1+normErrFraction)) 
            sfHists[process+"NormDown"].Add(hists[mcProcess][var], -(1/(1+normErrFraction))) 
        else: 
            if debugLevel > 0: print "Process",mcProcess,"has associated scale factors.  No further uncertainty will be applied to its normalization."
            sfHists[process+"NormUp"].Add(hists[mcProcess][var], -1) 
            sfHists[process+"NormDown"].Add(hists[mcProcess][var], -1) 

    #divide data/MC
    sfHists[process].Divide(hists[process][var])
    sfHists[process+"NormUp"].Divide(hists[process][var])
    sfHists[process+"NormDown"].Divide(hists[process][var])

    #zero any negative scale factors
    if isinstance(var, basestring) or len(var) == 1: #1D
        for bx in range(1, sfHists[process].GetNbinsX()+1):
            sfHists[process].SetBinContent(bx,max(0., sfHists[process].GetBinContent(bx)))
            sfHists[process+"NormUp"].SetBinContent(bx,max(0., sfHists[process+"NormUp"].GetBinContent(bx)))
            sfHists[process+"NormDown"].SetBinContent(bx,max(0., sfHists[process+"NormDown"].GetBinContent(bx)))
    elif len(var) == 2: #2D
        for bx in range(1, sfHists[process].GetNbinsX()+1):
            for by in range(1, sfHists[process].GetNbinsY()+1):
                sfHists[process].SetBinContent(bx,by,max(0., sfHists[process].GetBinContent(bx,by)))
                sfHists[process+"NormUp"].SetBinContent(bx,by,max(0., sfHists[process+"NormUp"].GetBinContent(bx,by)))
                sfHists[process+"NormDown"].SetBinContent(bx,by,max(0., sfHists[process+"NormDown"].GetBinContent(bx,by)))
    elif len(var) == 3: #3D
        for bx in range(1, sfHists[process].GetNbinsX()+1):
            for by in range(1, sfHists[process].GetNbinsY()+1):
                for bz in range(1, sfHists[process].GetNbinsZ()+1):
                    sfHists[process].SetBinContent(bx,by,bz,max(0., sfHists[process].GetBinContent(bx,by,bz)))
                    sfHists[process+"NormUp"].SetBinContent(bx,by,bz,max(0., sfHists[process+"NormUp"].GetBinContent(bx,by,bz)))
                    sfHists[process+"NormDown"].SetBinContent(bx,by,bz,max(0., sfHists[process+"NormDown"].GetBinContent(bx,by,bz)))

    #suppress scale factors consistent with 1
    if signifThreshold > 0:
        print "Ignoring scale factors compatible with 1.0 (",signifThreshold,"sigma significance )"
        if isinstance(var, basestring) or len(var) == 1: #1D
            for bx in range(1, sfHists[process].GetNbinsX()+1):
                if vetoLepCorrHist.GetBinError(bx) == 0: continue
                nsigma = abs(sfHists[process].GetBinContent(bx)-1.0)/sfHists[process].GetBinError(bx)
                if nsigma < signifThreshold: 
                    sfHists[process].SetBinContent(bx,1.0)
                    sfHists[process].SetBinError(bx,0.0)
        elif len(var) == 2: #2D
            for bx in range(1, sfHists[process].GetNbinsX()+1):
                for by in range(1, sfHists[process].GetNbinsY()+1):
                    nsigma = abs(sfHists[process].GetBinContent(bx,by)-1.0)/sfHists[process].GetBinError(bx,by)
                    if vetoLepCorrHist.GetBinError(bx,by) == 0: continue
                    if nsigma < signifThreshold: 
                        sfHists[process].SetBinContent(bx,by,1.0)
                        sfHists[process].SetBinError(bx,by,0.0)
        elif len(var) == 3: #3D
            for bx in range(1, sfHists[process].GetNbinsX()+1):
                for by in range(1, sfHists[process].GetNbinsY()+1):
                    for bz in range(1, sfHists[process].GetNbinsZ()+1):
                        if vetoLepCorrHist.GetBinError(bx,by,bz) == 0: continue
                        nsigma = abs(sfHists[process].GetBinContent(bx,by,bz)-1.0)/sfHists[process].GetBinError(bx,by,bz)
                        if nsigma < signifThreshold: 
                            sfHists[process].SetBinContent(bx,by,bz,1.0)
                            sfHists[process].SetBinError(bx,by,bz,0.0)

    if debugLevel > 0:
        print "Scale factor histograms after adding",process,":"
        print sfHists

    #plot scale factors in 2D (not yet implemented for 1 or 3 dimensions)
    c = rt.TCanvas("c"+process+"SFs", "c", 800, 600)
    if not isinstance(var, basestring) and len(var) == 2: 
        plotting.draw2DHist(c, sfHists[process], xtitle=var[0], ytitle=var[1], zmin=0.3, zmax=1.8, printstr=process+"ScaleFactors", lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", commentstr=process+" Data/MC Scale Factors", drawErrs=True, logz=False, numDigits=2, printdir=printdir)

        if printTable:
            xbinLowEdges = []
            xbinUpEdges = []
            ybinLowEdges = []
            ybinUpEdges = []
            sysUncerts = {mcProcess:[] for mcProcess in bgProcesses}
            statUncerts = []
            sfs = []
            #for each bin, get values for all table columns
            for bx in range(1, sfHists[process].GetNbinsX()+1):
                for by in range(1, sfHists[process].GetNbinsY()+1):
                    xbinLowEdges.append('%.0f' % (sfHists[process].GetXaxis().GetBinLowEdge(bx)))
                    xbinUpEdges.append('%.0f' % (sfHists[process].GetXaxis().GetBinUpEdge(bx)))
                    ybinLowEdges.append(str(sfHists[process].GetYaxis().GetBinLowEdge(by)))
                    ybinUpEdges.append(str(sfHists[process].GetYaxis().GetBinUpEdge(by)))
                    scaleFactor = sfHists[process].GetBinContent(bx,by)
                    sfs.append('%.3f' % (scaleFactor))
                    if scaleFactor > 0:
                        statUncerts.append('%.1f\\%%' % (100*(sfHists[process].GetBinError(bx,by)/scaleFactor)))
                    else: 
                        statUncerts.append('--')
                    for mcProcess in bgProcesses: 
                        dataYield = hists[dataName][var].GetBinContent(bx,by)
                        if dataYield > 0 and scaleFactor > 0:
                            sysUncerts[mcProcess].append('%.1f\\%%' % (100*abs(hists[mcProcess][var].GetBinContent(bx,by)*normErrFraction*1.0/dataYield/scaleFactor)))
                        else: 
                            sysUncerts[mcProcess].append('--')
            xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
            yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
            headers=[var[0], var[1], "Scale factor", "Stat.\\ unc."]
            cols = [xRanges, yRanges, sfs, statUncerts]
            for mcProcess in bgProcesses: 
                headers.extend(["Unc.\\ from "+mcProcess])
                cols.extend([sysUncerts[mcProcess]])
            plotting.table_basic(headers, cols, caption="Scale factors for "+process+" background", printstr="scaleFactorTable"+process, printdir=printdir)

def makeVetoLeptonCorrectionHist(hists={}, var=("MR","Rsq"), dataName="Data", lumiData=0, signifThreshold=0., debugLevel=0, regionName="Veto Lepton", printdir="."):
    """Compare data and MC in veto lepton control region.  Makes ratio histogram and returns it.  
    Arguments are similar to those for appendScaleFactors()
    """
    regionNameReduced = regionName.replace(' ','')

    #warn if the needed histograms are not found
    if dataName not in hists:
        print "Error in makeVetoLeptonCorrectionHist: target data histogram (",dataName,") was not found!"
        return
    if var not in hists[dataName]:
        print "Error in makeVetoLeptonCorrectionHist: could not find ",var," in hists[",dataName,"]!"

    #make total MC histogram
    bgProcesses = [mcProcess for mcProcess in hists if mcProcess != dataName and mcProcess != "Fit"] #get list of non-data, non-fit samples
    mcTotal = hists[dataName][var].Clone()
    mcTotal.Reset()
    for p in bgProcesses:
        if var not in hists[p]:
            print "Error in makeVetoLeptonCorrectionHist: histogram for",p,"not found!"
            return
        mcTotal.Add(hists[p][var])

    #subtract MC-data
    print "Making correction histogram for veto lepton control region"
    vetoLepCorrHist = mcTotal.Clone(regionNameReduced+"Correction")
    vetoLepCorrHist.SetDirectory(0)
    vetoLepCorrHist.Add(hists[dataName][var],-1.0)

    #suppress corrections consistent with 0
    if signifThreshold > 0:
        print "Ignoring corrections compatible with 0 at",signifThreshold,"sigma significance"
        if isinstance(var, basestring) or len(var) == 1: #1D
            for bx in range(1, vetoLepCorrHist.GetNbinsX()+1):
                if vetoLepCorrHist.GetBinError(bx) == 0: continue
                nsigma = abs(vetoLepCorrHist.GetBinContent(bx))/vetoLepCorrHist.GetBinError(bx)
                if nsigma < signifThreshold: 
                    vetoLepCorrHist.SetBinContent(bx,0.0)
                    vetoLepCorrHist.SetBinError(bx,0.0)
        elif len(var) == 2: #2D
            for bx in range(1, vetoLepCorrHist.GetNbinsX()+1):
                for by in range(1, vetoLepCorrHist.GetNbinsY()+1):
                    if vetoLepCorrHist.GetBinError(bx,by) == 0: continue
                    nsigma = abs(vetoLepCorrHist.GetBinContent(bx,by))/vetoLepCorrHist.GetBinError(bx,by)
                    if nsigma < signifThreshold: 
                        vetoLepCorrHist.SetBinContent(bx,by,0.0)
                        vetoLepCorrHist.SetBinError(bx,by,0.0)
        elif len(var) == 3: #3D
            for bx in range(1, vetoLepCorrHist.GetNbinsX()+1):
                for by in range(1, vetoLepCorrHist.GetNbinsY()+1):
                    for bz in range(1, vetoLepCorrHist.GetNbinsZ()+1):
                        if vetoLepCorrHist.GetBinError(bx,by,bz) == 0: continue
                        nsigma = abs(vetoLepCorrHist.GetBinContent(bx,by,bz))/vetoLepCorrHist.GetBinError(bx,by,bz)
                        if nsigma < signifThreshold: 
                            vetoLepCorrHist.SetBinContent(bx,by,bz,0.0)
                            vetoLepCorrHist.SetBinError(bx,by,bz,0.0)

    #plot correction factors in 2D (not yet implemented for 1 or 3 dimensions)
    c = rt.TCanvas("cVetoLeptonCorrs", "c", 800, 600)
    if not isinstance(var, basestring) and len(var) == 2: 
        plotting.draw2DHist(c, vetoLepCorrHist, xtitle=var[0], ytitle=var[1], zmin=-200, zmax=200, printstr=regionNameReduced+"Correction", lumistr=('%.1f' % (lumiData*1.0/1000))+" fb^{-1}", commentstr="MC-Data, "+regionName+" Control Region", drawErrs=True, logz=False, numDigits=2, printdir=printdir)

        xbinLowEdges = []
        xbinUpEdges = []
        ybinLowEdges = []
        ybinUpEdges = []
        statUncerts = []
        sfs = []
        #for each bin, get values for all table columns
        for bx in range(1, vetoLepCorrHist.GetNbinsX()+1):
            for by in range(1, vetoLepCorrHist.GetNbinsY()+1):
                xbinLowEdges.append('%.0f' % (vetoLepCorrHist.GetXaxis().GetBinLowEdge(bx)))
                xbinUpEdges.append('%.0f' % (vetoLepCorrHist.GetXaxis().GetBinUpEdge(bx)))
                ybinLowEdges.append(str(vetoLepCorrHist.GetYaxis().GetBinLowEdge(by)))
                ybinUpEdges.append(str(vetoLepCorrHist.GetYaxis().GetBinUpEdge(by)))
                corrFactor = vetoLepCorrHist.GetBinContent(bx,by)
                sfs.append('%.3f' % (corrFactor))
                if corrFactor != 0:
                    statUncerts.append('%.1f\\%%' % (100*(vetoLepCorrHist.GetBinError(bx,by)/corrFactor)))
                else: 
                    statUncerts.append('--')
                for mcProcess in bgProcesses: 
                    dataYield = hists[dataName][var].GetBinContent(bx,by)
        xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
        yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
        headers=[var[0], var[1], "MC-Data", "Stat.\\ unc."]
        cols = [xRanges, yRanges, sfs, statUncerts]
        plotting.table_basic(headers, cols, caption="Difference between MC and data yields in veto lepton control region", printstr="corrFactorTable"+regionNameReduced, printdir=printdir)

    return vetoLepCorrHist
