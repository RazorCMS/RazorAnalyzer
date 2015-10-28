## Machinery for inclusive razor analysis

import sys,os
from array import array
from collections import namedtuple
import ROOT as rt

#local imports
import macro
from razorAnalysis import *
from razorWeights import *

###########################################
### BASIC HISTOGRAM FILLING/PLOTTING MACRO
###########################################

def makeControlSampleHists(regionName="TTJetsSingleLepton", filenames={}, samples=[], cutsMC="", cutsData="", bins={}, logX=True, lumiMC=1, lumiData=3000, weightHists={}, sfHists={}, treeName="ControlSampleEvent",dataName="Data", weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], shapeErrors=[], miscErrors=[], fitToyFiles=None, boxName="", btags=-1, blindBins=None, debugLevel=0):
    titles = {
        "MR": "M_{R} (GeV)", 
        "Rsq": "R^{2}",
        "mll": "m_{ll} (GeV)",
        }
    #setup files and trees
    inputs = filenames
    files = {name:rt.TFile.Open(inputs[name]) for name in inputs} #get input files
    for name in inputs: 
        assert files[name] #check input files
        if debugLevel > 0: print "Opened file",inputs[name]
    trees = macro.makeTreeDict(files, treeName, debugLevel)

    #add noise filters to cuts if bits are present in tree
    cutsData = appendNoiseFilters(cutsData, trees[dataName]) 

    #define histograms to fill
    hists,shapeHists = macro.makeHistograms(regionName, inputs, samples, bins, titles, shapeErrors, dataName)
    listOfVars = hists.itervalues().next().keys() #list of the variable names
    
    #fill histograms
    print("\nData:")
    macro.loopTree(trees[dataName], weightF=weight_data, cuts=cutsData, hists=hists[dataName], weightHists=weightHists, weightOpts=[], debugLevel=debugLevel) 

    print("\nMC:")
    macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:hists[name] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, debugLevel=debugLevel) 

    #get up/down histogram variations
    for shape in shapeErrors:
        print "\n"+shape,"Up:"
        macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:shapeHists[name][shape+"Up"] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, errorOpt=shape+"Up", boxName=boxName, debugLevel=debugLevel)
        print "\n"+shape,"Down:"
        macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:shapeHists[name][shape+"Down"] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, weightOpts=weightOpts, errorOpt=shape+"Down", boxName=boxName, debugLevel=debugLevel)

    #propagate up/down systematics to central histograms
    for name in samples:
        for var in bins:
            for shape in shapeErrors:
                if debugLevel > 0: print "Adding",shape,"uncertainty in quadrature with",name,"errors for",var
                #loop over histogram bins
                for bx in range(1, hists[name][var].GetNbinsX()+1):
                    #use difference between Up and Down histograms as uncertainty
                    sysErr = abs(shapeHists[name][shape+'Up'][var].GetBinContent(bx) - shapeHists[name][shape+'Down'][var].GetBinContent(bx))
                    #add in quadrature with existing error
                    oldErr = hists[name][var].GetBinError(bx)
                    hists[name][var].SetBinError(bx, (oldErr**2 + sysErr**2)**(0.5))
                    if debugLevel > 0: print shape,": Error on bin ",bx,"increases from",oldErr,"to",hists[name][var].GetBinError(bx),"after adding",sysErr,"in quadrature"
            for source in miscErrors:
                if source.lower() == "mt" and var == "MR":
                    applyMTUncertainty1D(hists[name][var], process=name+"_"+boxName, debugLevel=debugLevel)
        #2D case
        if "MR" in bins and "Rsq" in bins:
            for shape in shapeErrors:
                if debugLevel > 0: print "Adding",shape,"uncertainty in quadrature with",name,"errors for razor histogram"
                #loop over histogram bins
                for bx in range(1, hists[name][("MR","Rsq")].GetNbinsX()+1):
                    for by in range(1, hists[name][("MR","Rsq")].GetNbinsY()+1):
                    #use difference between Up and Down histograms as uncertainty
                        sysErr = abs(shapeHists[name][shape+'Up'][("MR","Rsq")].GetBinContent(bx,by) - shapeHists[name][shape+'Down'][("MR","Rsq")].GetBinContent(bx,by))
                        #add in quadrature with existing error
                        oldErr = hists[name][("MR","Rsq")].GetBinError(bx,by)
                        hists[name][("MR","Rsq")].SetBinError(bx,by, (oldErr**2 + sysErr**2)**(0.5))
                        if debugLevel > 0: print shape,": Error on bin (",bx,by,") increases from",oldErr,"to",hists[name][("MR","Rsq")].GetBinError(bx,by),"after adding",sysErr,"in quadrature"
            for source in miscErrors:
                if source.lower() == "mt":
                    applyMTUncertainty2D(hists[name][("MR","Rsq")], process=name+"_"+boxName, debugLevel=debugLevel)

    c = rt.TCanvas(regionName+"c", regionName+"c", 800, 600)

    #optionally compare with fit
    if fitToyFiles and "MR" in bins and "Rsq" in bins:
        print "Comparing MC prediction with fit"
        #make histogram for fit result
        hists["Fit"] = {}
        for v in ["MR","Rsq",("MR","Rsq")]:
            hists["Fit"][v] = hists[dataName][v].Clone(hists[dataName][v].GetName()+"Fit")
        hists["Fit"][v].Reset()

        #load fit information, including toys
        toyFile = rt.TFile(fitToyFiles[boxName])
        assert toyFile
        if debugLevel > 0: print "Opened file",fitToyFiles[boxName],"to get Bayesian toy results"
        toyTree = toyFile.Get("myTree")
        assert toyTree
        if debugLevel > 0: print "Got tree myTree"

        #get uncertainties on fit prediction
        #using code imported from Javier's plotting script
        from PlotFit import getBinSumDicts, getBestFitRms, getErrors1D
        z = [0.,1.,2.,3.,4.] #b-tag binning
        if btags < 0:
            zmin = 0
            zmax = len(z)-1
        else:
            zmin = btags+1
            zmax = btags+1

        #make dummy options tuple
        Opt = namedtuple("Opt", "printErrors")
        options = Opt(False)

        #store best fit and uncertainty in each bin
        #1D
        binSumDict = getBinSumDicts('x', 0,len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
        print binSumDict
        for i,sumName in binSumDict.iteritems():
            nObs = hists[dataName]["MR"].GetBinContent(i)
            bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options,"")
            hists["Fit"]["MR"].SetBinContent(i,bestFit)
            hists["Fit"]["MR"].SetBinError(i,rms)
        binSumDict = getBinSumDicts('y', 0, len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
        for i,sumName in binSumDict.iteritems():
            nObs = hists[dataName]["Rsq"].GetBinContent(i)
            bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options,"")
            hists["Fit"]["Rsq"].SetBinContent(i,bestFit)
            hists["Fit"]["Rsq"].SetBinError(i,rms)
        #2D
        binSumDict = getBinSumDicts('yx', 0, len(bins['MR'])-1,0,len(bins['Rsq'])-1,zmin,zmax,bins['MR'],bins['Rsq'],z)
        for (i,j),sumName in binSumDict.iteritems():
            nObs = hists[dataName][("MR","Rsq")].GetBinContent(i,j)
            bestFit, rms, pvalue, nsigma, c = getBestFitRms(toyTree,sumName,nObs,c,options=options,plotName="")
            hists["Fit"][("MR","Rsq")].SetBinContent(i,j,bestFit)
            hists["Fit"][("MR","Rsq")].SetBinError(i,j,rms)

    #print histograms
    rt.SetOwnership(c, False)
    macro.basicPrint(hists, mcNames=samples, varList=listOfVars, c=c, printName=regionName, logx=logX, dataName=dataName, ymin=0.1, lumistr=str(lumiData)+" pb^{-1}", boxName=boxName, btags=btags, blindBins=blindBins)

    #close files and return
    for f in files: files[f].Close()
    return hists

#######################################
### MAKE SCALE FACTORS FROM HISTOGRAMS
#######################################

def appendScaleFactors(process="TTJets", hists={}, sfHists={}, var=("MR","Rsq"), dataName="Data", normErrFraction=0.2, printTable=True, lumiData=0, debugLevel=0):
    """Subtract backgrounds and make the data/MC histogram for the given process.
    Also makes up/down histograms corresponding to 20% shifts in the background normalization"""
    if debugLevel > 0: 
        print "Scale factor histograms so far:"
        print sfHists
    #warn if this scale factor histogram already exists
    if process in sfHists:
        print "Warning in appendScaleFactors: ",process," scale factor histogram already exists!  Will overwrite..."
    #get the target MC histogram
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
    print "Making scale factor histogram for",process
    sfHists[process] = hists[dataName][var].Clone(process+"ScaleFactors")
    sfHists[process].SetDirectory(0)
    #make systematic error histograms
    sfHists[process+"NormUp"] = hists[dataName][var].Clone(process+"ScaleFactorsUp")
    sfHists[process+"NormUp"].SetDirectory(0)
    sfHists[process+"NormDown"] = hists[dataName][var].Clone(process+"ScaleFactorsDown")
    sfHists[process+"NormDown"].SetDirectory(0)

    #subtract backgrounds in data
    bgProcesses = [mcProcess for mcProcess in hists if mcProcess != process and mcProcess != dataName and mcProcess != "Fit"]
    for mcProcess in bgProcesses:
        if var not in hists[mcProcess]:
            print "Error in appendScaleFactors: could not find",var," in hists[",mcProcess,"]!"
            return
        if debugLevel > 0: print "Subtracting",mcProcess,"from",dataName,"distribution"
        sfHists[process].Add(hists[mcProcess][var], -1) 
        if mcProcess not in sfHists: #background normalization uncertainty affects only processes with no scale factors
            if debugLevel > 0: print "Process",mcProcess,"has no associated scale factors.  Its normalization will be given a 20% uncertainty"
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
    for bx in range(1, sfHists[process].GetNbinsX()+1):
        for by in range(1, sfHists[process].GetNbinsY()+1):
            sfHists[process].SetBinContent(bx,by,max(0., sfHists[process].GetBinContent(bx,by)))
            sfHists[process+"NormUp"].SetBinContent(bx,by,max(0., sfHists[process+"NormUp"].GetBinContent(bx,by)))
            sfHists[process+"NormDown"].SetBinContent(bx,by,max(0., sfHists[process+"NormDown"].GetBinContent(bx,by)))

    if debugLevel > 0:
        print "Scale factor histograms after adding",process,":"
        print sfHists

    #plot scale factors in 2D
    c = rt.TCanvas("c"+process+"SFs", "c", 800, 600)
    macro.draw2DHist(c, sfHists[process], xtitle=var[0], ytitle=var[1], zmin=0.3, zmax=1.8, printstr=process+"ScaleFactors", lumistr=str(lumiData)+" pb^{-1}", commentstr=process+" Data/MC Scale Factors", drawErrs=True)

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
        macro.table_basic(headers, cols, caption="Scale factors for "+process+" background", printstr="scaleFactorTable"+process)
