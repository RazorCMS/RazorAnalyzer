import os
import ROOT as rt
import copy
from array import *

#local imports
from PlotFit import setFFColors
from RunCombine import exec_me
from razorWeights import applyMTUncertainty1D, applyMTUncertainty2D

def makeGrayGraphs(hNS):    
    if hNS is None or hNS == 0: return
    fGrayGraphs = []
    col1 = rt.gROOT.GetColor(rt.kGray+1)
    col1.SetAlpha(0.3)
    for iBinX in range(1,hNS.GetNbinsX()+1):
        for iBinY in range(1,hNS.GetNbinsY()+1):
            if hNS.GetBinContent(iBinX,iBinY)!= -999: continue
            xBinLow = hNS.GetXaxis().GetBinLowEdge(iBinX)
            xBinHigh = xBinLow+hNS.GetXaxis().GetBinWidth(iBinX)
            yBinLow = hNS.GetYaxis().GetBinLowEdge(iBinY)
            yBinHigh = yBinLow+hNS.GetYaxis().GetBinWidth(iBinY)
            fGray = rt.TGraph(5)
            fGray.SetPoint(0,xBinLow,yBinLow)
            fGray.SetPoint(1,xBinLow,yBinHigh)
            fGray.SetPoint(2,xBinHigh,yBinHigh)
            fGray.SetPoint(3,xBinHigh,yBinLow)
            fGray.SetPoint(4,xBinLow,yBinLow)
            fGray.SetFillColor(rt.kGray+1)
            fGrayGraphs.append(fGray)
    return fGrayGraphs

def blindHistograms(histList, blindBins):
    """blindBins should be a list of ordered pairs corresponding to bin coordinates"""
    for hist in histList:
        for (x,y) in blindBins:
            hist.SetBinContent(x,y,-999)

def setupHistograms(regionName, inputs, samples, bins, titles, shapeErrors, dataName):
    hists = {name:{} for name in inputs}
    shapeHists = {name:{} for name in inputs}
    if inputs is None: return
    for name in inputs:
        #1D histograms
        for var in bins:
            if var in titles: title=titles[var]
            else: title = var
            hists[name][var] = rt.TH1F(regionName+var+name, title, len(bins[var])-1, array('d',bins[var]))
            #add up/down histograms for each systematic uncertainty
            if samples is not None and name in samples:
                for shape in shapeErrors:
                    if shape+"Down" not in shapeHists[name]: shapeHists[name][shape+"Down"] = {}
                    if shape+"Up" not in shapeHists[name]: shapeHists[name][shape+"Up"] = {}
                    shapeHists[name][shape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Down")
                    shapeHists[name][shape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Up")
        #2D MR-Rsq histogram
        if "MR" in bins and "Rsq" in bins:
            hists[name][("MR","Rsq")] = rt.TH2F(regionName+"MRRsq"+name, "R^{2} vs M_{R}", len(bins["MR"])-1, array('d',bins["MR"]), len(bins["Rsq"])-1, array('d',bins["Rsq"]))
            if samples is not None and name in samples: 
                for shape in shapeErrors:
                    shapeHists[name][shape+"Down"][("MR","Rsq")] = hists[name][("MR","Rsq")].Clone(hists[name][("MR","Rsq")].GetName()+shape+"Down")
                    shapeHists[name][shape+"Up"][("MR","Rsq")] = hists[name][("MR","Rsq")].Clone(hists[name][("MR","Rsq")].GetName()+shape+"Up")
        for var in hists[name]: 
            hists[name][var].Sumw2()
            hists[name][var].SetDirectory(0)
            if samples is not None and name in samples:
                for shape in shapeErrors:
                    shapeHists[name][shape+"Down"][var].Sumw2()
                    shapeHists[name][shape+"Up"][var].Sumw2()
                    shapeHists[name][shape+"Down"][var].SetDirectory(0)
                    shapeHists[name][shape+"Up"][var].SetDirectory(0)
    #deal correctly with Poisson errors on data
    if dataName in hists:
        for var in hists[dataName]:
            hists[dataName][var].SetBinErrorOption(rt.TH1.kPoisson)
    return hists,shapeHists

def propagateShapeSystematics(hists, samples, bins, shapeHists, shapeErrors, miscErrors=[], boxName="", debugLevel=0):
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

def basicPrint(histDict, mcNames, varList, c, printName="Hist", dataName="Data", logx=False, ymin=0.1, lumistr="40 pb^{-1}", boxName=None, btags=None, comment=True, blindBins=None, nsigmaFitData=None, nsigmaFitMC=None, doDensity=False, printdir="."):
    """Make stacked plots of quantities of interest, with data overlaid"""
    #format MC histograms
    for name in mcNames: 
        for var in histDict[name]: setHistColor(histDict[name][var], name)
    titles = {name:name for name in mcNames}

    #get data histograms
    if dataName in histDict: dataHists = histDict[dataName]
    else: dataHists = None

    #make correct comment string
    if comment:
        if boxName is None or boxName == "":
            commentstr = printName
        else:
            commentstr = boxName+ " Box"
            if btags is not None and btags >= 0:
                commentstr += ", "+str(btags)+" B-tag"
    else:
        commentstr = ""

    legend=None
    plotFit = ("Fit" in histDict)
    for i,var in enumerate(varList): 
        #for MR and Rsq, make 2D plots
        if var == ('MR','Rsq'):
            if len(mcNames) > 0:
                mcPrediction = histDict[mcNames[0]][var].Clone(histDict[mcNames[0]][var].GetName()+"mcPrediction")
                mcPrediction.Reset()
                for name in mcNames: 
                    mcPrediction.Add(histDict[name][var])
            else:
                mcPrediction = 0
            #copy data and fit histograms
            if dataHists is not None: obsData = dataHists[var].Clone("obsData")
            else: obsData = 0
            if not plotFit: 
                fitPrediction = 0
            else: 
                fitPrediction = histDict["Fit"][var]
            #blind signal region if necessary
            if blindBins is not None and dataHists is not None:
                blindHistograms([obsData], blindBins)
                if nsigmaFitData is not None:
                    blindHistograms([nsigmaFitData], blindBins)
            #make plots
            plot_basic_2D(c, mc=mcPrediction, data=obsData, fit=fitPrediction, xtitle='MR', ytitle='Rsq', printstr='Razor_'+printName, lumistr=lumistr, commentstr=commentstr, saveroot=True, nsigmaFitData=nsigmaFitData, nsigmaFitMC=nsigmaFitMC, printdir=printdir)
            #print prediction in each bin
            if obsData != 0:
                print "Results for razor data histogram:"
                for bx in range(1, obsData.GetNbinsX()+1):
                    for by in range(1,obsData.GetNbinsY()+1):
                        print bx,by,obsData.GetBinContent(bx,by),"+/-",obsData.GetBinError(bx,by)
                print "\n"
            if mcPrediction != 0:
                print "Results for razor MC histogram:"
                for bx in range(1, mcPrediction.GetNbinsX()+1):
                    for by in range(1,mcPrediction.GetNbinsY()+1):
                        print bx,by,mcPrediction.GetBinContent(bx,by),"+/-",mcPrediction.GetBinError(bx,by)
                print "\n"
            if fitPrediction != 0:
                print "Results for razor fit histogram:"
                for bx in range(1, fitPrediction.GetNbinsX()+1):
                    for by in range(1,fitPrediction.GetNbinsY()+1):
                        print bx,by,fitPrediction.GetBinContent(bx,by),"+/-",fitPrediction.GetBinError(bx,by)
                print "\n"
        #for other variables make 1D plots
        if not isinstance(var, basestring): continue #only consider strings
        varHists = {name:histDict[name][var].Clone(histDict[name][var].GetName()+"Clone") for name in mcNames}
        if doDensity: #divide each histogram bin by its width
            for name in varHists:
                for bx in range(1,varHists[name].GetNbinsX()+1):
                    varHists[name].SetBinContent(bx, varHists[name].GetBinContent(bx)*1.0/varHists[name].GetXaxis().GetBinWidth(bx))
                    varHists[name].SetBinError(bx, varHists[name].GetBinError(bx)*1.0/varHists[name].GetXaxis().GetBinWidth(bx))
        stack = makeStack(varHists, mcNames, var)
        if len(mcNames) == 0: stack = None #plotting function can't handle an empty stack
        if dataHists is not None:
            obsData = dataHists[var].Clone(dataHists[var].GetName()+"obsData")
            if doDensity:
                for bx in range(1,obsData.GetNbinsX()+1):
                    obsData.SetBinContent(bx, obsData.GetBinContent(bx)*1.0/obsData.GetXaxis().GetBinWidth(bx))
                    obsData.SetBinError(bx, obsData.GetBinError(bx)*1.0/obsData.GetXaxis().GetBinWidth(bx))
        else:
            obsData = None
        if not plotFit:
            fitPrediction = None
        else:
            fitPrediction = histDict["Fit"][var].Clone(histDict["Fit"][var].GetName()+"FitPrediction")
            if doDensity:
                for bx in range(1,fitPrediction.GetNbinsX()+1):
                    fitPrediction.SetBinContent(bx, fitPrediction.GetBinContent(bx)*1.0/fitPrediction.GetXaxis().GetBinWidth(bx))
                    fitPrediction.SetBinError(bx, fitPrediction.GetBinError(bx)*1.0/fitPrediction.GetXaxis().GetBinWidth(bx))
        if not legend:
            legend = makeLegend(varHists, titles, reversed(mcNames))
            if blindBins is None and obsData is not None: legend.AddEntry(obsData, dataName)
            if plotFit and fitPrediction is not None: legend.AddEntry(fitPrediction, "Fit")
        if doDensity:
            ytitle = "Events / Bin Width"
            ymin = None
        else:
            ytitle = "Events"
        if blindBins is None:
            plot_basic(c, mc=stack, data=obsData, fit=fitPrediction, leg=legend, xtitle=var, ytitle=ytitle, printstr=var+"_"+printName, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, printdir=printdir)
        else:
            plot_basic(c, mc=stack, data=None, fit=fitPrediction, leg=legend, xtitle=var, ytitle=ytitle, printstr=var+"_"+printName, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, printdir=printdir)

def transformVarsInString(string, varNames, suffix):
    outstring = copy.copy(string)
    for var in varNames:
        outstring = outstring.replace(var, var+suffix) 
    return outstring

def transformVarString(event, string, errorOpt, process="", debugLevel=0):
    jetvars = ["MR","Rsq","nBTaggedJets","dPhiRazor","leadingJetPt","subleadingJetPt","nSelectedJets","nJets80","box"] #quantities susceptible to jet uncertainties
    outstring = string
    if errorOpt == "jesUp":
        outstring = transformVarsInString(string, jetvars, "_JESUp")
    elif errorOpt == "jesDown":
        oustring = transformVarsInString(string, jetvars, "_JESDown")
    elif errorOpt == "jerUp":
        outstring = transformVarsInString(string, jetvars, "_JERUp")
    elif errorOpt == "jerDown":
        outstring = transformVarsInString(string, jetvars, "_JERDown")

    if debugLevel > 1:
        if outstring != string: print "For option",errorOpt,"Replacing string '",string,"' with '",outstring,"'"
    return outstring

def getAdditionalCuts(tree, errorOpt, process, debugLevel=0):
    """Implement here any event-by-event cuts that can't be handled with TTree::Draw"""
    pass

def basicFill(tree, hists={}, weight=1.0, sysErrSquaredHists={}, sysErr=0.0, errorOpt=None, additionalCuts=None, debugLevel=0):
    """Fills each histogram with the corresponding variable in the tree.
    'hists' should be a dictionary of histograms, with keys being the variable names to fill.
    Ex: hists['MR'] should be the histogram you want to fill with MR values.
    A key that is a tuple of variables (ex: ('MR','Rsq')) should be paired with a multidimensional histogram.
    In this case, the given variables will be used to fill the histogram."""
    #make additional cuts 
    if additionalCuts is not None:
        additionalCutsBool = eval(additionalCuts)
        if debugLevel > 1: 
            print "Additional cuts:",additionalCuts,"\nDecision:",additionalCutsBool
        if not additionalCutsBool: return
    #fill each variable
    for varName, hist in hists.iteritems(): 
        if isinstance(varName, basestring): #if varName is a string
            #transform variable name
            if errorOpt is not None: varName = transformVarString(tree, varName, errorOpt, debugLevel=debugLevel)
            if debugLevel > 1: print "Filling",varName,"=",getattr(tree,varName),"with weight",weight
            hist.Fill(getattr(tree, varName), weight)
            if varName in sysErrSquaredHists: #for propagating systematic errors on the variables
                sysErrSquared = weight*weight*sysErr*sysErr
                sysErrSquaredHist[varName].Fill(getattr(tree, varName), sysErrSquared)
                if debugLevel > 1: print "Sys. Error =",sysErr,";Filling (w*sysErr)^2 histogram with",sysErrSquared
        else: #treat it as a tuple of variables that should be filled
            #transform each variable
            if errorOpt is not None: varName = tuple([transformVarString(tree, v, errorOpt, debugLevel=debugLevel) for v in varName])
            toFill = [getattr(tree, v) for v in varName]+[weight]
            if debugLevel > 1: print "Filling",varName,":",toFill
            hist.Fill(*toFill)
            if varName in sysErrSquaredHists:
                sysErrSquared = weight*weight*sysErr*sysErr
                toFillErr = [getattr(tree, v) for v in varName]+[sysErrSquared]
                sysErrSquaredHists[varName].Fill(*toFillErr)
                if debugLevel > 1: print "Sys. Error =",sysErr,";Filling (w*sysErr)^2 histogram with",sysErrSquared

def makeTreeDict(fileDict, treeName, debugLevel=0):
    """gets a tree called treeName from each file in fileDict, and returns a dict of trees"""
    trees = {}
    for name in fileDict:
        if debugLevel > 0: print("Loading tree "+treeName+" for process "+name)
        trees[name] = fileDict[name].Get(treeName)
        if debugLevel > 0: print("Got tree containing "+str(trees[name].GetEntries())+" entries")
        assert trees[name]
    if debugLevel > 0: 
        print("Trees loaded: ") 
        print trees
    return trees

def getScaleFactorAndError(tree, sfHist, sfVars=("MR","Rsq"), debugLevel=0):
    #get variables
    var = [getattr(tree, v) for v in sfVars]
    #constrain variables to be within the bounds of the histogram
    var[0] = min(var[0], sfHist.GetXaxis().GetXmax()*0.999)
    var[0] = max(var[0], sfHist.GetXaxis().GetXmin()*1.001)
    if len(var) > 1:
        var[1] = min(var[1], sfHist.GetYaxis().GetXmax()*0.999)
        var[1] = max(var[1], sfHist.GetYaxis().GetXmin()*1.001)
    if len(var) > 2:
        var[2] = min(var[2], sfHist.GetZaxis().GetXmax()*0.999)
        var[2] = max(var[2], sfHist.GetZaxis().GetXmin()*1.001)
    scaleFactor = sfHist.GetBinContent(sfHist.FindFixBin(*var))
    scaleFactorErr = sfHist.GetBinError(sfHist.FindFixBin(*var))
    if debugLevel > 1: print "Applying scale factor: ",scaleFactor
    return (scaleFactor, scaleFactorErr)

def addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel=0):
    """For each histogram in hists, look for the corresponding histogram in sysErrSquaredHists.
    Treats the values of sysErrSquaredHists as sums of (weight*error)^2 in each bin, and adds these errors in quadrature with the existing bin errors in hists"""
    for name in hists:
        if name in sysErrSquaredHists:
            if debugLevel > 0: print "Including systematic errors on ",name
            for bx in range(1, hists[name].GetNbinsX()+1):
                for by in range(1, hists[name].GetNbinsY()+1):
                    squaredError = sysErrSquaredHists[name].GetBinContent(bx,by)
                    oldErr = hists[name].GetBinError(bx,by)
                    hists[name].SetBinError(bx,by,(oldErr*oldErr + squaredError)**(0.5))
                    if debugLevel > 0: print name,": Error on bin (",bx,by,") increases from",oldErr,"to",hists[name].GetBinError(bx,by),"after adding",(squaredError**(0.5)),"in quadrature"

def loopTree(tree, weightF, cuts="", hists={}, weightHists={}, sfHist=None, scale=1.0, fillF=basicFill, sfVars=("MR","Rsq"), sysVars=("MR", "Rsq"), weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], errorOpt=None, process="", debugLevel=0):
    """Loop over a single tree and fill histograms.
    Returns the sum of the weights of selected events."""
    if debugLevel > 0: print ("Looping tree "+tree.GetName())
    if debugLevel > 0: print ("Cuts: "+cuts)
    #transform cuts 
    if errorOpt is not None:
        if debugLevel > 0: print "Error option is:",errorOpt
        cuts = transformVarString(tree, cuts, errorOpt, process=process, debugLevel=debugLevel+1)
    #get list of entries passing the cuts
    tree.Draw('>>elist', cuts, 'entrylist')
    elist = rt.gDirectory.Get('elist')
    print "Total entries passing cuts:",elist.GetN()
    #create histograms for systematics
    sysErrSquaredHists = {}
    for name in hists: 
        if name == sysVars:
            sysErrSquaredHists[name] = hists[name].Clone(hists[name].GetName()+"ERRORS")
            sysErrSquaredHists[name].Reset()
            if debugLevel > 0: print "Created temp histogram",sysErrSquaredHists[name].GetName(),"to hold",name,"systematic errors"
    count = 0
    sumweight = 0.0
    while True:
        #load the next entry
        entry = elist.Next()
        if entry == -1: break
        if count > 0 and count % 100000 == 0: print "Processing entry",count
        elif debugLevel > 0 and count % 10000 == 0: print "Processing entry",count
        elif debugLevel > 1: print "Processing entry",count
        tree.GetEntry(entry)
        w = weightF(tree, weightHists, scale, weightOpts, errorOpt, debugLevel=debugLevel)
        additionalCuts = getAdditionalCuts(tree, errorOpt, process, debugLevel=debugLevel) #additional cuts according to systematic
        err = 0.0
        if sfHist is not None: 
            sf, err = getScaleFactorAndError(tree, sfHist, sfVars, debugLevel)
            w *= sf
        fillF(tree, hists, w, sysErrSquaredHists, err, errorOpt, additionalCuts, debugLevel)
        sumweight += w
        count += 1
    #propagate systematics to each histogram
    addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel)
    print "Sum of weights for this sample:",sumweight
    return sumweight

def loopTrees(treeDict, weightF, cuts="", hists={}, weightHists={}, sfHists={}, scale=1.0, weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], errorOpt=None, fillF=basicFill, sfVars=("MR","Rsq"), sysVars=("MR","Rsq"), boxName="NONE", debugLevel=0):
    """calls loopTree on each tree in the dictionary.  
    Here hists should be a dict of dicts, with hists[name] the collection of histograms to fill using treeDict[name]"""
    sumweights=0.0
    for name in treeDict: 
        if name not in hists: continue
        print("Filling histograms for tree "+name)
        sfHistToUse = None
        if name in sfHists: 
            print("Using scale factors from histogram "+sfHists[name].GetName())
            sfHistToUse = sfHists[name]
        sumweights += loopTree(treeDict[name], weightF, cuts, hists[name], weightHists, sfHistToUse, scale, fillF, sfVars, sysVars, weightOpts, errorOpt, process=name+"_"+boxName, debugLevel=debugLevel)
    print "Sum of event weights for all processes:",sumweights

def makeStack(hists, ordering, title="Stack"):
    """Takes a dict of histograms and an ordered list of names, and returns a THStack containing the histograms stacked in the desired order"""
    stack = rt.THStack("thstack"+title.replace(" ",""), title)
    rt.SetOwnership(stack, False)
    for name in ordering: 
        if name in hists:
            stack.Add(hists[name])
        else: 
            print("Warning in makeStack: histogram "+name+" not found in histogram dictionary!")
    return stack

def makeLegend(hists, titles, ordering, x1=0.6, y1=0.6, x2=0.9, y2=0.9):
    """Takes a dict of histograms, a dict of histogram titles, and an ordered list of names, and returns a legend with the histograms in the desired order"""
    leg = rt.TLegend(x1, y1, x2, y2)
    rt.SetOwnership(leg, False)
    for name in ordering: 
        leg.AddEntry(hists[name], titles[name])
    return leg

colors = {"WJets":rt.kRed+1, "DYJets":rt.kBlue+1, "TTJets":rt.kGreen+2, "ZInv":rt.kCyan+1, "QCD":rt.kOrange+3, "SingleTop":rt.kOrange-3, "VV":rt.kViolet+3, "TTV":rt.kGreen-7, "DYJetsLow":rt.kBlue+1, "Other":rt.kViolet+3}
def setHistColor(hist, name):
    """Sets histogram color"""
    if name in colors: hist.SetFillColor(colors[name])
    else: print "Warning in macro.py: histogram fill color not set for",name

def plot_basic(c, mc=0, data=0, fit=0, leg=0, xtitle="", ytitle="Number of events", ymin=None, ymax=None, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", commentstr="", ratiomin=0.5, ratiomax=1.5, pad2Opt="Ratio", fitColor=rt.kBlue, mcErrColor=rt.kRed, customPad2Hist=None, saveroot=False, savepdf=False, savepng=True, printdir='.'):
    """Plotting macro with options for data, MC, and fit histograms.  Creates data/MC ratio if able."""
    #setup
    c.Clear()
    c.cd()
    if (data and mc) or (data and fit) or (mc and fit): pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.4, 1, 1)
    else: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.1, 1, 1)
    pad1.SetBottomMargin(0)
    pad1.SetLogx(logx)
    pad1.SetLogy(logy)
    pad1.Draw()
    pad1.cd()

    #draw MC
    if mc:
        #make total MC histogram
        histList = mc.GetHists()
        mcTotal = histList.First().Clone()
        mcTotal.SetStats(0)
        mcTotal.Reset()
        for h in histList:
            mcTotal.Add(h)
        mcTotal.SetFillColor(mcErrColor-10)
        mcTotal.SetTitle("")
        mcTotal.GetYaxis().SetTitle(ytitle)
        mcTotal.GetYaxis().SetLabelSize(0.03)
        if logy: mcTotal.GetXaxis().SetMoreLogLabels()
        if not data: mcTotal.GetXaxis().SetTitle(xtitle)
        if data: mcTotal.GetYaxis().SetTitleOffset(0.45)
        else: mcTotal.GetYaxis().SetTitleOffset(0.50)
        mcTotal.GetYaxis().SetTitleSize(0.05)
        if ymin is not None: mcTotal.SetMinimum(ymin)
        if ymax is not None: mcTotal.SetMaximum(ymax)
        if data and data.GetMaximum() > mcTotal.GetMaximum() and ymax is None: 
            mcTotal.SetMaximum(data.GetMaximum())
        mcTotal.Draw("e2")
        mc.Draw("histsame")
    #draw fit
    if fit:
        fit.SetStats(0)
        fit.SetMarkerStyle(21)
        fit.SetLineWidth(2)
        fit.SetLineColor(fitColor+2)
        fit.SetMarkerColor(fitColor)
        fit.GetYaxis().SetTitle(ytitle)
        fit.SetTitle("")
        if ymin is not None: fit.SetMinimum(ymin)
        if ymax is not None: fit.SetMaximum(ymax)
        fit.Draw("pesame")
    #draw data
    if data:
        data.SetStats(0)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1)
        data.GetYaxis().SetTitle(ytitle)
        data.SetTitle("")
        if ymin is not None: data.SetMinimum(ymin)
        if ymax is not None: data.SetMaximum(ymax)
        data.SetBinErrorOption(rt.TH1.kPoisson)
        data.Draw("pe1same")
    pad1.Modified()
    rt.gPad.Update()

    #add legend and LaTeX 
    leg.Draw()
    t1 = rt.TLatex(0.1,0.94, "CMS Preliminary")
    t2 = rt.TLatex(0.55,0.94, "#sqrt{s}=13 TeV"+((lumistr != "")*(", L = "+lumistr)))
    t1.SetNDC()
    t2.SetNDC()
    t1.SetTextSize(0.06)
    t2.SetTextSize(0.06)
    t1.Draw()
    t2.Draw()
    if commentstr != "":
        t3 = rt.TLatex(0.30, 0.84, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.04)
        t3.Draw()

    #lower pad plot
    lowerPadHist = None

    #set up custom histogram if provided
    if customPad2Hist is not None:
        customPad2Hist.SetTitle("")
        customPad2Hist.GetXaxis().SetTitle(xtitle)
        customPad2Hist.SetMinimum(ratiomin)
        customPad2Hist.SetMaximum(ratiomax)
        customPad2Hist.SetStats(0)
        customPad2Hist.GetXaxis().SetLabelSize(0.1)
        customPad2Hist.GetYaxis().SetLabelSize(0.08)
        customPad2Hist.GetYaxis().SetTitleOffset(0.35)
        customPad2Hist.GetXaxis().SetTitleOffset(1.00)
        customPad2Hist.GetYaxis().SetTitleSize(0.08)
        customPad2Hist.GetXaxis().SetTitleSize(0.08)
        customPad2Hist.SetTitle("")
        if logx: customPad2Hist.GetXaxis().SetMoreLogLabels()

    #make ratio data/MC
    if data and mc:
        if customPad2Hist is not None:
            lowerPadHist = customPad2Hist.Clone("lowerPadHist")
            lowerPadHist.GetYaxis().SetTitle("Nsigma")
            lowerPadHist.SetMarkerStyle(20)
            lowerPadHist.SetMarkerSize(1)
        elif pad2Opt.lower() == "ratio":
            lowerPadHist = make1DRatioHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("Data / MC")
        elif pad2Opt.lower() == "nsigma" or pad2Opt.lower() == "pulls" or pad2Opt.lower() == "ff":
            lowerPadHist = make1DPullHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - MC)/#sigma")
        elif pad2Opt.lower() == "perc" or pad2Opt.lower() == "percent":
            lowerPadHist = make1DPercentDiffHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - MC)/MC")
    #make ratio data/fit
    elif data and fit:
        if customPad2Hist is not None:
            lowerPadHist = customPad2Hist.Clone("lowerPadHist")
            lowerPadHist.GetYaxis().SetTitle("Nsigma")
            lowerPadHist.SetMarkerStyle(20)
            lowerPadHist.SetMarkerSize(1)
        elif pad2Opt.lower() == "ratio":
            lowerPadHist = make1DRatioHistogram(data, fit, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("Data / Fit")
        elif pad2Opt.lower() == "nsigma" or pad2Opt.lower() == "pulls" or pad2Opt.lower() == "ff":
            lowerPadHist = make1DPullHistogram(data, fit, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - Fit)/#sigma")
        elif pad2Opt.lower() == "perc" or pad2Opt.lower() == "percent":
            lowerPadHist = make1DPercentDiffHistogram(data, fit, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Data - Fit)/Fit")
    #make ratio mc/fit
    elif mc and fit:
        if customPad2Hist is not None:
            lowerPadHist = customPad2Hist.Clone("lowerPadHist")
            lowerPadHist.GetYaxis().SetTitle("Nsigma")
            lowerPadHist.SetMarkerStyle(fit.GetMarkerStyle())
            lowerPadHist.SetMarkerSize(fit.GetMarkerSize())
        elif pad2Opt.lower() == "ratio":
            lowerPadHist = make1DRatioHistogram(fit, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("Fit / MC")
        elif pad2Opt.lower() == "nsigma" or pad2Opt.lower() == "pulls" or pad2Opt.lower() == "ff":
            lowerPadHist = make1DPullHistogram(fit, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Fit - MC)/#sigma")
        elif pad2Opt.lower() == "perc" or pad2Opt.lower() == "percent":
            lowerPadHist = make1DPercentDiffHistogram(fit, mcTotal, xtitle, ratiomin, ratiomax, logx)
            lowerPadHist.GetYaxis().SetTitle("(Fit - MC)/Fit")
    c.cd()
    pad2 = rt.TPad(printstr+"pad2",printstr+"pad2",0,0.0,1,0.4)
    pad2.SetTopMargin(0)
    pad2.SetTopMargin(0.008)
    pad2.SetBottomMargin(0.25)
    pad2.SetGridy()
    pad2.SetLogx(logx)
    if lowerPadHist is not None:
        pad2.Draw()
        pad2.cd()
        lowerPadHist.Draw("pe")
        pad2.Modified()
        rt.gPad.Update()

    #save
    if savepng: c.Print(printdir+'/'+printstr+".png")
    if savepdf: c.Print(printdir+'/'+printstr+".pdf")
    if saveroot: c.Print(printdir+'/'+printstr+".root")

def draw2DHist(c, hist, xtitle="", ytitle="", ztitle="", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, drawErrs=False, palette=53, grayGraphs=None, saveroot=False, savepdf=False, savepng=True, numDigits=1, printdir='.'):
    """Draw a single 2D histogram and print to file"""
    if palette == "FF":
        setFFColors(hist, -5.1, 5.1)
    else:
        rt.gStyle.SetPalette(palette)
    c.Clear()
    c.cd()
    c.SetLogx(logx)
    c.SetLogy(logy)
    c.SetLogz(logz)
    c.Draw()
    c.cd()
    hist.SetTitle("")
    hist.GetXaxis().SetTitle(xtitle)
    hist.GetYaxis().SetTitle(ytitle)
    #hist.GetZaxis().SetTitle(ztitle)
    hist.GetZaxis().SetTitle("") #until we can get the z-axis to display correctly
    hist.GetYaxis().SetLabelSize(0.03)
    hist.GetYaxis().SetTitleOffset(0.50)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.SetStats(0)
    if zmin is not None: hist.SetMinimum(zmin)
    elif hist.GetMinimum() < -10: hist.SetMinimum(0.1) #avoid drawing -999's etc
    if zmax is not None: hist.SetMaximum(zmax)
    hist.Draw("colz")
    if grayGraphs is not None: 
        for g in grayGraphs: g.Draw("f")
    if dotext:
        rt.gStyle.SetPaintTextFormat('4.%df' % numDigits)
        hist.SetMarkerSize(2.0)
        if not drawErrs: hist.Draw('textsame')
        else: 
            hist.SetMarkerSize(1.0)
            hist.Draw('textesame')
    #add LaTeX 
    t1 = rt.TLatex(0.1,0.94, "CMS Preliminary")
    t2 = rt.TLatex(0.55,0.94, "#sqrt{s}=13 TeV"+((lumistr != "")*(", L = "+lumistr)))
    t1.SetNDC()
    t2.SetNDC()
    t1.SetTextSize(0.05)
    t2.SetTextSize(0.05)
    t1.Draw()
    t2.Draw()
    if commentstr != "":
        t3 = rt.TLatex(0.30, 0.84, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.04)
        t3.Draw()
    #save
    if savepng: c.Print(printdir+'/'+printstr+".png")
    if savepdf: c.Print(printdir+'/'+printstr+".pdf")
    if saveroot: c.Print(printdir+'/'+printstr+".root")

def make1DRatioHistogram(num, denom, xtitle="", ratiomin=0.25, ratiomax=2.0, logx=False, forPad2=True):
    ratio = num.Clone()
    ratio.Divide(denom)
    ratio.SetTitle("")
    ratio.GetXaxis().SetTitle(xtitle)
    ratio.SetMinimum(ratiomin)
    ratio.SetMaximum(ratiomax)
    ratio.SetStats(0)
    if forPad2:
        ratio.GetXaxis().SetLabelSize(0.1)
        ratio.GetYaxis().SetLabelSize(0.08)
        ratio.GetYaxis().SetTitleOffset(0.35)
        ratio.GetXaxis().SetTitleOffset(1.00)
        ratio.GetYaxis().SetTitleSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.08)
        ratio.SetTitle("")
        if logx: ratio.GetXaxis().SetMoreLogLabels()
    return ratio

def make1DPullHistogram(h1, h2, xtitle="", ymin=-5.0, ymax=5.0, logx=False, forPad2=True, suppress=True, suppressLevel=0.1):
    """Makes (h1 - h2)/sigma histogram, where sigma is the error on the difference"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"Pulls")
    ret.Add(h2, -1)
    for bx in range(1, h1.GetNbinsX()+1):
        content = ret.GetBinContent(bx)
        err1 = h1.GetBinError(bx)
        err2 = h2.GetBinError(bx)
        err = (err1*err1 + err2*err2)**(0.5)
        if err > 0:
            ret.SetBinContent(bx,content*1.0/err)
            ret.SetBinError(bx, 0.0)
        else:
            ret.SetBinContent(bx,-9999)
        if suppress:
            if h1.GetBinContent(bx) < suppressLevel and h2.GetBinContent(bx) < suppressLevel:
                ret.SetBinContent(bx,-9999)
    ret.SetMinimum(ymin)
    ret.SetMaximum(ymax)
    if forPad2:
        ret.GetXaxis().SetLabelSize(0.1)
        ret.GetYaxis().SetLabelSize(0.08)
        ret.GetYaxis().SetTitleOffset(0.35)
        ret.GetXaxis().SetTitleOffset(1.00)
        ret.GetYaxis().SetTitleSize(0.08)
        ret.GetXaxis().SetTitleSize(0.08)
        ret.SetTitle("")
        if logx: ret.GetXaxis().SetMoreLogLabels()

    return ret

def make1DPercentDiffHistogram(h1, h2, xtitle="", ymin=-1.0, ymax=1.0, logx=False, forPad2=True):
    """Makes (h1 - h2)/h2 histogram"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"PercDiff")
    ret.Add(h2, -1)
    ret.Divide(h2)
    ret.SetMinimum(ymin)
    ret.SetMaximum(ymax)
    if forPad2:
        ret.GetXaxis().SetLabelSize(0.1)
        ret.GetYaxis().SetLabelSize(0.08)
        ret.GetYaxis().SetTitleOffset(0.35)
        ret.GetXaxis().SetTitleOffset(1.00)
        ret.GetYaxis().SetTitleSize(0.08)
        ret.GetXaxis().SetTitleSize(0.08)
        ret.SetTitle("")
        if logx: ret.GetXaxis().SetMoreLogLabels()
    return ret

def make2DPullHistogram(h1, h2, suppress=True):
    """Makes (h1 - h2)/sigma histogram, where sigma is the error on the difference"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"Pulls")
    ret.Add(h2, -1)
    for bx in range(1, h1.GetNbinsX()+1):
        for by in range(1, h1.GetNbinsY()+1):
            content = ret.GetBinContent(bx,by)
            err1 = h1.GetBinError(bx,by)
            err2 = h2.GetBinError(bx,by)
            err = (err1*err1 + err2*err2)**(0.5)
            if err > 0:
                ret.SetBinContent(bx,by,content*1.0/err)
            else:
                ret.SetBinContent(bx,by,0)
            #suppress bins that have < 1 entry in both histograms
            if suppress and h1.GetBinContent(bx,by) < 1 and h2.GetBinContent(bx,by) < 1:
                ret.SetBinContent(bx,by,-9999)
    return ret

def make2DPercentDiffHistogram(h1, h2, suppress=True):
    """Makes (h1 - h2)/h2 histogram"""
    ret = h1.Clone(h1.GetName()+h2.GetName()+"PercentDiff")
    ret.Add(h2, -1)
    ret.Divide(h2)
    if suppress:
        for bx in range(1, h1.GetNbinsX()+1):
            for by in range(1, h1.GetNbinsY()+1):
                #suppress bins that have < 1 entry in both histograms
                if h1.GetBinContent(bx,by) < 1 and h2.GetBinContent(bx,by) < 1:
                    ret.SetBinContent(bx,by,-9999)
    return ret

def make2DRelativeUncertaintyHistogram(h, suppress=True, suppressLevel=10.0):
    ret = h.Clone(h.GetName()+"RelUnc")
    for bx in range(1, h.GetNbinsX()+1):
        for by in range(1, h.GetNbinsY()+1):
            if h.GetBinContent(bx,by) != 0:
                ret.SetBinContent(bx,by,h.GetBinError(bx,by)*1.0/h.GetBinContent(bx,by))
            else:
                ret.SetBinContent(bx,by,-9999)
            if suppress and ret.GetBinContent(bx,by) > suppressLevel:
                ret.SetBinContent(bx,by,-9999)
    return ret

def unroll2DHistograms(hists):
    out = [] 
    for hist in hists:
        if hist is None or hist == 0: 
            out.append(None)
            continue
        numbins = hist.GetNbinsX()*hist.GetNbinsY()
        outHist = rt.TH1F(hist.GetName()+"Unroll", hist.GetTitle(), numbins, 0, numbins)
        ibin = 0
        for bx in range(1, hist.GetNbinsX()+1):
            for by in range(1, hist.GetNbinsY()+1):
                ibin += 1
                outHist.SetBinContent(ibin, hist.GetBinContent(bx,by))
                outHist.SetBinError(ibin, hist.GetBinError(bx,by))
        out.append(outHist)
    return out

def plot_basic_2D(c, mc=0, data=0, fit=0, xtitle="", ytitle="", ztitle="Number of events", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, saveroot=False, savepdf=False, savepng=True, nsigmaFitData=None, nsigmaFitMC=None, printdir="."):
    """Plotting macro for data, MC, and/or fit yields.  Creates french flag plots comparing data/MC/fit if able."""
    #make a gray square for each -999 bin
    grayGraphs = [makeGrayGraphs(hist) for hist in [mc,fit,data]]
    #unroll 2D hists to plot in 1D
    unrolled = unroll2DHistograms([mc, data, fit])
    if data: unrolled[1].SetBinErrorOption(rt.TH1.kPoisson) #get correct error bars on data
    #synchronize z-axes
    if zmin is None:
        smallestZ = 0.1
        if mc: smallestZ = min(smallestZ, max(0.1,mc.GetMinimum()/1.2))
        if data: smallestZ = min(smallestZ, max(0.1,data.GetMinimum()/1.2))
        if fit: smallestZ = min(smallestZ, max(0.1,fit.GetMinimum()/1.2))
        if mc: mc.SetMinimum(smallestZ)
        if data: data.SetMinimum(smallestZ)
        if fit: fit.SetMinimum(smallestZ)
    if zmax is None:
        largestZ = -float('inf')
        if mc: largestZ = max(largestZ, mc.GetMaximum()*1.2)
        if data: largestZ = max(largestZ, data.GetMaximum()*1.2)
        if fit: largestZ = max(largestZ, fit.GetMaximum()*1.2)
        if mc: mc.SetMaximum(largestZ)
        if data: data.SetMaximum(largestZ)
        if fit: fit.SetMaximum(largestZ)
    #draw each histogram and all relevant combinations
    if mc:
        unrolled[0].SetLineColor(rt.kBlack)
        unrolled[0].SetFillColor(38)
        mcStack = makeStack({"MC":unrolled[0]}, ["MC"], "MC")
        draw2DHist(c, mc, xtitle, ytitle, ztitle, zmin, zmax, printstr+'MC', lumistr=lumistr, commentstr=commentstr+", MC prediction", dotext=dotext, grayGraphs=grayGraphs[0], drawErrs=True, saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
        if data: 
            #do (data - mc)/unc
            mcPulls = make2DPullHistogram(data,mc)
            draw2DHist(c, mcPulls, xtitle, ytitle, ztitle, None, None, printstr+'MCPulls', lumistr=lumistr, commentstr=commentstr+", (Data - MC)/#sigma", dotext=dotext, palette="FF", logz=False, grayGraphs=grayGraphs[2], saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
            #do (data - mc)/mc
            mcPerc = make2DPercentDiffHistogram(data,mc)
            draw2DHist(c, mcPerc, xtitle, ytitle, ztitle, -1.5, 1.5, printstr+'MCPercentDiff', lumistr=lumistr, commentstr=commentstr+", (Data - MC)/MC", palette="FF", dotext=dotext, logz=False, grayGraphs=grayGraphs[2], saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
            #unroll and compare
            blindMC = unrolled[0].Clone("blindMC")
            for bx in range(1,blindMC.GetNbinsX()+1):
                if unrolled[1].GetBinContent(bx) < 0:
                    blindMC.SetBinContent(bx,-999)
            blindStack = makeStack({"MC":blindMC}, ["MC"], "MC")
            legDataMC = rt.TLegend(0.7, 0.7, 0.9, 0.9)
            rt.SetOwnership(legDataMC, False)
            legDataMC.AddEntry(blindMC, "MC Prediction")
            legDataMC.AddEntry(unrolled[1], "Data")
            plot_basic(c, blindStack, unrolled[1], None, legDataMC, xtitle="Bin", ymin=0.1, printstr=printstr+"UnrolledDataMC", lumistr=lumistr, commentstr=commentstr, ratiomin=-5.0,ratiomax=5.0, pad2Opt="ff", saveroot=True, printdir=printdir)
        if fit: 
            #do (fit - mc)/unc
            if nsigmaFitMC is None:
                #do (MC - fit)/unc
                mcFitPulls = make2DPullHistogram(fit,mc)
                note="(MC - Fit)/#sigma"
                printnote="MCFitPulls"
            else:
                #make nsigma plot
                note="Nsigmas"
                printnote="MCFitNSigma"
                mcFitPulls = nsigmaFitMC
            draw2DHist(c, mcFitPulls, xtitle, ytitle, ztitle, None, None, printstr+printnote, lumistr=lumistr, commentstr=commentstr+note, palette="FF", logz=False, dotext=dotext, grayGraphs=grayGraphs[0], saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
            #do (fit - mc)/mc
            mcFitPerc = make2DPercentDiffHistogram(fit,mc)
            draw2DHist(c, mcFitPerc, xtitle, ytitle, ztitle, -1.5, 1.5, printstr+'MCFitPercentDiff', lumistr=lumistr, commentstr=commentstr+", (Fit - MC)/MC", palette="FF", logz=False, dotext=dotext, grayGraphs=grayGraphs[0], saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
            #unroll and compare

            nsigmaUnrolledFitMC = None
            if nsigmaFitMC is not None:
                nsigmaUnrolledFitMC = unroll2DHistograms([nsigmaFitMC])[0]
                for bx in range(1, nsigmaUnrolledFitMC.GetNbinsX()+1):
                    nsigmaUnrolledFitMC.SetBinError(bx,0.0)

            legMCFit = rt.TLegend(0.7, 0.7, 0.9, 0.9)
            rt.SetOwnership(legMCFit, False)
            legMCFit.AddEntry(unrolled[0], "MC Prediction")
            legMCFit.AddEntry(unrolled[2], "Fit Prediction")
            plot_basic(c, mcStack, None, unrolled[2], legMCFit, xtitle="Bin", ymin=0.1, printstr=printstr+"UnrolledMCFit", lumistr=lumistr, commentstr=commentstr, ratiomin=-5., ratiomax=5.0, pad2Opt="ff", saveroot=True, mcErrColor=rt.kRed, customPad2Hist=nsigmaUnrolledFitMC, printdir=printdir)
    if data:
        draw2DHist(c, data, xtitle, ytitle, ztitle, zmin=max(0.1,zmin), printstr=printstr+'Data', lumistr=lumistr, commentstr=commentstr+", Data", dotext=dotext, grayGraphs=grayGraphs[2], saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
        if fit: 
            if nsigmaFitData is None:
                #do (data - fit)/unc
                dataFitPulls = make2DPullHistogram(data,fit)
                note=", (Data - Fit)/#sigma"
                printnote="DataFitPulls"
            else:
                #make nsigma plot
                note=", Nsigmas"
                printnote="DataFitNSigma"
                dataFitPulls = nsigmaFitData
            draw2DHist(c, dataFitPulls, xtitle, ytitle, ztitle, None, None, printstr+printnote, lumistr=lumistr, commentstr=commentstr+note, dotext=dotext, palette="FF", logz=False, grayGraphs=grayGraphs[2], saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
            #do (data - fit)/fit
            dataFitPerc = make2DPercentDiffHistogram(data,fit)
            draw2DHist(c, dataFitPerc, xtitle, ytitle, ztitle, -1.5, 1.5, printstr+'DataFitPercentDiff', lumistr=lumistr, commentstr=commentstr+", (Data - Fit)/Fit", palette="FF", dotext=dotext, logz=False, grayGraphs=grayGraphs[2], saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
            #unroll and compare
            blindFit = unrolled[2].Clone("blindFit")

            nsigmaUnrolled = None
            if nsigmaFitData is not None:
                nsigmaUnrolled = unroll2DHistograms([nsigmaFitData])[0]
                for bx in range(1,nsigmaUnrolled.GetNbinsX()+1):
                    nsigmaUnrolled.SetBinError(bx,0.0)

            for bx in range(1,blindFit.GetNbinsX()+1):
                if unrolled[1].GetBinContent(bx) < 0:
                    blindFit.SetBinContent(bx, -999)

            legDataFit = rt.TLegend(0.7, 0.7, 0.9, 0.9)
            rt.SetOwnership(legDataFit, False)
            legDataFit.AddEntry(blindFit, "Fit Prediction")
            legDataFit.AddEntry(unrolled[1], "Data")

            plot_basic(c, None, unrolled[1], blindFit, legDataFit, xtitle="Bin", ymin=0.1, printstr=printstr+"UnrolledDataFit", lumistr=lumistr, commentstr=commentstr, ratiomin=-5., ratiomax=5.0, pad2Opt="ff", fitColor=rt.kGreen, saveroot=True, customPad2Hist=nsigmaUnrolled, printdir=printdir)
    if fit:
        draw2DHist(c, fit, xtitle, ytitle, ztitle, zmin, zmax, printstr+'Fit', lumistr=lumistr, commentstr=commentstr+", Fit prediction", grayGraphs=grayGraphs[1], dotext=dotext, drawErrs=True, saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
        #make relative uncertainty histogram
        relUncFitHist = make2DRelativeUncertaintyHistogram(fit)
        draw2DHist(c, relUncFitHist, xtitle, ytitle, ztitle, 0.0, 2.0, printstr+'FitRelUnc', lumistr=lumistr, commentstr=commentstr+", Fit relative uncertainty", grayGraphs=grayGraphs[1], dotext=dotext, drawErrs=False, saveroot=saveroot, savepdf=savepdf, savepng=savepng, logz=False, printdir=printdir)

def makeStackAndPlot(canvas, mcHists={}, dataHist=None, dataName="Data", mcOrdering=[], titles=[], mcTitle="Stack", xtitle="", ytitle="Number of events", printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", saveroot=False, savepdf=False, savepng=True, ymin=None, ymax=None):
    #make stack
    stack = makeStack(mcHists, mcOrdering, mcTitle)
    #make legend
    hists = copy.copy(mcHists)
    hists[dataName] = dataHist
    ordering = copy.copy(mcOrdering)
    ordering.append(dataName)
    leg = makeLegend(hists, titles, ordering)
    #plot
    plot_basic(canvas, stack, dataHist, leg=leg, xtitle=xtitle, ytitle=ytitle, printstr=printstr, logx=logx, logy=logy, lumistr=lumistr, saveroot=saveroot, savepdf=savepdf, savepng=savepng, ymin=ymin)

def table_basic(headers=[], cols=[], caption="", printstr='table', landscape=False):
    #check for input
    if len(cols) == 0:
        print "table_basic: no columns provided.  doing nothing."
        return
    #check that all columns have the same length
    for col in cols:
        if len(col) != len(cols[0]):
            print "Error in table_basic: columns do not have equal lengths!"
            return
    #check that there is a header for each column
    if len(headers) != len(cols):
        print "Error in table_basic: number of headers does not equal number of columns!"
        return

    with open(printstr+'.tex', 'w') as f:
        #f.write('\\newgeometry{margin=0.2cm}\n')
        if landscape: f.write('\\begin{landscape}\n')
        f.write('\\begin{center}\n\\footnotesize\n\\begin{longtable}{|'+('|'.join(['c' for c in cols]))+'|}\n')
        f.write('\\caption{'+caption+'}\n\\endhead\n\\hline\n')
        f.write(' & '.join(headers)+' \\\\\n\\hline\n')
        for row in range(len(cols[0])):
            f.write((' & '.join([col[row] for col in cols]))+' \\\\\n\\hline\n')
        f.write('\\end{longtable}\n\\end{center}\n')
        if landscape: f.write('\\end{landscape}\n')
        #f.write('\\restoregeometry\n')
        print "Created LaTeX scale factor table",(printstr+".tex")
