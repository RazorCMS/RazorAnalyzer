import os
import ROOT as rt
import copy
import math
from array import *

#local imports
from PlotFit import setFFColors
from RunCombine import exec_me
from razorWeights import applyMTUncertainty1D, applyMTUncertainty2D
from plotting import *

def blindHistograms(histList, blindBins):
    """blindBins should be a list of ordered pairs corresponding to bin coordinates.
    Sets blinded bin contents to -999."""
    for hist in histList:
        for (x,y) in blindBins:
            hist.SetBinContent(x,y,-999)

def setupHistograms(regionName, inputs, samples, bins, titles, shapeErrors, dataName):
    """Creates dictionary of histograms with specified binning.

    regionName: used only in histogram names and titles
    inputs: dictionary of process:filename pairs
    samples: list of processes for which histograms should be made
    bins: dictionary.  key = name of observable, value = list of bins for that observable
    titles: dictionary.  key = name of observable, value = title for histograms of that quantity 
    shapeErrors: list of shape uncertainties to apply
    """

    hists = {name:{} for name in inputs}
    shapeHists = {name:{} for name in inputs}
    if inputs is None: return
    for name in inputs:
        for var in bins:
            if isinstance(var, basestring): 
                #1D histograms
                if var in titles: title=titles[var]
                else: title = var
                hists[name][var] = rt.TH1F(regionName+var+name, title+';'+title+';', len(bins[var])-1, array('d',bins[var]))
                #add up/down histograms for each systematic uncertainty
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if shape+"Down" not in shapeHists[name]: shapeHists[name][shape+"Down"] = {}
                        if shape+"Up" not in shapeHists[name]: shapeHists[name][shape+"Up"] = {}
                        shapeHists[name][shape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Down")
                        shapeHists[name][shape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Up")
            elif len(var) == 2: 
                #2D histograms
                title = [titles[v] if v in titles else v for v in var]
                hists[name][var] = rt.TH2F(regionName+var[0]+var[1]+name, ';'+title[0]+';'+title[1], len(bins[var[0]])-1, array('d',bins[var[0]]), len(bins[var[1]])-1, array('d',bins[var[1]]))
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if shape+"Down" not in shapeHists[name]: shapeHists[name][shape+"Down"] = {}
                        if shape+"Up" not in shapeHists[name]: shapeHists[name][shape+"Up"] = {}
                        shapeHists[name][shape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Down")
                        shapeHists[name][shape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Up")
            elif len(var) == 3:
                #3D histograms
                title = [titles[v] if v in titles else v for v in var]
                hists[name][var] = rt.TH3F(regionName+var[0]+var[1]+var[2]+name, ';'+title[0]+';'+title[1]+';'+title[2], len(bins[var[0]])-1, array('d',bins[var[0]]), len(bins[var[1]])-1, array('d',bins[var[1]]), len(bins[var[2]])-1, array('d',bins[var[2]]))
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if shape+"Down" not in shapeHists[name]: shapeHists[name][shape+"Down"] = {}
                        if shape+"Up" not in shapeHists[name]: shapeHists[name][shape+"Up"] = {}
                        shapeHists[name][shape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Down")
                        shapeHists[name][shape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+shape+"Up")
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
    """For each bin of the central histogram, add the appropriate uncertainties in quadrature with the statistical uncertainty.
    List of arguments is similar to razorMacros.makeControlSampleHists
    NOTE: not set up for 2D or 3D histograms yet!"""

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
                    if debugLevel > 0 and sysErr > 0: print shape,": Error on bin ",bx,"increases from",oldErr,"to",hists[name][var].GetBinError(bx),"after adding",sysErr,"in quadrature"
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

def basicPrint(histDict, mcNames, varList, c, printName="Hist", dataName="Data", logx=False, ymin=0.1, lumistr="40 pb^{-1}", boxName=None, btags=None, comment=True, blindBins=None, nsigmaFitData=None, nsigmaFitMC=None, doDensity=False, printdir=".", special="", vartitles={}):
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
            commentstr = '#it{ razor '+boxName
            if btags is not None and btags >= 0:
                commentstr += " "+str(btags)+" b-tag"
            if 'sideband' in special:
                commentstr += " Sideband Fit"
            elif 'full' in special:
                commentstr += ' Full Fit'
            commentstr += '}'
    else:
        commentstr = ""

    legend=None
    plotFit = ("Fit" in histDict)
    for i,var in enumerate(varList): 
        #for MR and Rsq, make 2D plots
        if not isinstance(var, basestring) and len(var) == 2: #2D plots
            mcDict = None 
            if len(mcNames) > 0:
                mcDict = {} #for stacked unrolled plots
                mcPrediction = histDict[mcNames[0]][var].Clone(histDict[mcNames[0]][var].GetName()+"mcPrediction")
                mcPrediction.Reset()
                for name in mcNames: 
                    mcPrediction.Add(histDict[name][var])
                    mcDict[name] = histDict[name][var] #for stacked unrolled plots
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
            #get axis titles
            if var[0] in vartitles:
                xtitle = vartitles[var[0]]
            else:
                xtitle = var[0]
            if var[1] in vartitles:
                ytitle = vartitles[var[1]]
            else:
                ytitle = var[1]
            #make plots
            plot_basic_2D(c, mc=mcPrediction, data=obsData, fit=fitPrediction, xtitle=xtitle, ytitle=ytitle, printstr=var[0]+var[1]+printName, lumistr=lumistr, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, nsigmaFitData=nsigmaFitData, nsigmaFitMC=nsigmaFitMC, mcDict=mcDict, mcSamples=mcNames, printdir=printdir)
            #print prediction in each bin
            if obsData != 0:
                print "Results for data histogram:"
                for bx in range(1, obsData.GetNbinsX()+1):
                    for by in range(1,obsData.GetNbinsY()+1):
                        print bx,by,obsData.GetBinContent(bx,by),"+/-",obsData.GetBinError(bx,by)
                print "\n"
            if mcPrediction != 0:
                print "Results for MC histogram:"
                for bx in range(1, mcPrediction.GetNbinsX()+1):
                    for by in range(1,mcPrediction.GetNbinsY()+1):
                        print bx,by,mcPrediction.GetBinContent(bx,by),"+/-",mcPrediction.GetBinError(bx,by)
                print "\n"
            if fitPrediction != 0:
                print "Results for fit histogram:"
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
        #set x title
        if var in vartitles:
            xtitle = vartitles[var]
        else:
            xtitle = var
        if var in ['MR','Rsq','MR_NoW',"Rsq_NoW","MR_NoZ","Rsq_NoZ", "lep1.Pt()"]: logx = True
        else: logx = False
        if blindBins is None:
            plot_basic(c, mc=stack, data=obsData, fit=fitPrediction, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=var+"_"+printName, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, printdir=printdir)
        else:
            plot_basic(c, mc=stack, data=None, fit=fitPrediction, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=var+"_"+printName, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, printdir=printdir)

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

def basicFill(tree, hists={}, weight=1.0, sysErrSquaredHists={}, sysErr=0.0, errorOpt=None, additionalCuts=None, formulas={}, debugLevel=0):
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
            if varName in formulas: #if TTreeFormula
                varValue = formulas[varName].EvalInstance()
            else: #if tree variable
                varValue =  getattr(tree, varName)
            hist.Fill(varValue, weight)
            if debugLevel > 1: print "Filling",varName,"=",varValue,"with weight",weight
            if varName in sysErrSquaredHists: #for propagating systematic errors on the variables
                sysErrSquared = weight*weight*sysErr*sysErr
                sysErrSquaredHist[varName].Fill(varValue, sysErrSquared)
                if debugLevel > 1: print "Sys. Error =",sysErr,";Filling (w*sysErr)^2 histogram with",sysErrSquared
        else: #treat it as a tuple of variables that should be filled
            #transform each variable
            if errorOpt is not None: varName = tuple([transformVarString(tree, v, errorOpt, debugLevel=debugLevel) for v in varName])
            toFill = [formulas[v].EvalInstance() if v in formulas else getattr(tree, v) for v in varName]+[weight]
            if debugLevel > 1: print "Filling",varName,":",toFill
            hist.Fill(*toFill)
            if varName in sysErrSquaredHists:
                sysErrSquared = weight*weight*sysErr*sysErr
                toFillErr = [formulas[v].EvalInstance() if v in formulas else getattr(tree, v) for v in varName]+[sysErrSquared]
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

def getScaleFactorAndError(tree, sfHist, sfVars=("MR","Rsq"), formulas={}, debugLevel=0):
    if isinstance(sfVars, basestring):
        sfVars = (sfVars,) #cast into tuple 
    #get variables
    var = [formulas[v].EvalInstance() if v in formulas else getattr(tree, v) for v in sfVars]
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

    #add protection against unphysics scale factors
    #if scale factor is 0, then use scale factor of 1 and add 100% systematic uncertainty
    if scaleFactor == 0:
        scaleFactor = 1.0
        scaleFactorErr = math.sqrt( 1 + scaleFactorErr*scaleFactorErr )
    #if scale factor is > 2.0, then use scale factor of 1 and add difference between scale factor and 1.0 as a systematic uncertainty
    if scaleFactor > 2.0:
        scaleFactorErr = math.sqrt( (scaleFactor-1) + scaleFactorErr*scaleFactorErr )
        scaleFactor = 1.0        

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

def loopTree(tree, weightF, cuts="", hists={}, weightHists={}, sfHist=None, scale=1.0, fillF=basicFill, sfVars=("MR","Rsq"), sysVars=("MR", "Rsq"), weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], errorOpt=None, process="", auxSFs={}, auxSFHists={}, debugLevel=0):
    """Loop over a single tree and fill histograms.
    Returns the sum of the weights of selected events."""
    if debugLevel > 0: print ("Looping tree "+tree.GetName())
    if debugLevel > 0: print ("Cuts: "+cuts)
    #make TTreeFormulas for variables not found in the tree
    formulas = {}
    for var in hists:
        if isinstance(var, basestring) and not hasattr(tree, var): #if it's not in the tree
            formulas[var] = rt.TTreeFormula(var, var, tree)
            formulas[var].GetNdata()
    #make TTreeFormulas for auxSFs
    auxSFForms = {}
    for name,pair in auxSFs.iteritems(): #auxSFs pairs look like "ScaleFactor":("Var to reweight","Cuts")
        if debugLevel > 0: print "Making TTreeFormula for",name,"with formula",pair[1],"for reweighting",pair[0]
        auxSFForms[name] = rt.TTreeFormula(name+"Cuts", pair[1], tree)
        auxSFForms[name].GetNdata()
    #transform cuts 
    if errorOpt is not None:
        if debugLevel > 0: print "Error option is:",errorOpt
        cuts = transformVarString(tree, cuts, errorOpt, process=process, debugLevel=debugLevel+1)
    #get list of entries passing the cuts
    tree.Draw('>>elist', cuts, 'entrylist')
    elist = rt.gDirectory.Get('elist')
    print "Total entries passing cuts:",elist.GetN()
    #create histograms for scale factor systematics
    sysErrSquaredHists = {}
    for name in hists: 
        if name == sysVars:
            sysErrSquaredHists[name] = hists[name].Clone(hists[name].GetName()+"SFERRORS")
            sysErrSquaredHists[name].Reset()
            if debugLevel > 0: print "Created temp histogram",sysErrSquaredHists[name].GetName(),"to hold systematic errors from",name,"scale factors"
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
        additionalCuts = getAdditionalCuts(tree, errorOpt, process, debugLevel=debugLevel) #currently does nothing
        err = 0.0
        if sfHist is not None: 
            sf, err = getScaleFactorAndError(tree, sfHist, sfVars, formulas, debugLevel)
            w *= sf
        for name in auxSFs: #apply misc scale factors (e.g. veto lepton correction)
            if auxSFForms[name].EvalInstance(): #check if this event should be reweighted
                auxSF, auxErr = getScaleFactorAndError(tree, auxSFHists[name], sfVars=auxSFs[name][0], formulas=formulas, debugLevel=debugLevel)
                w *= auxSF
                err = (err*err + auxErr*auxErr)**(0.5)
        fillF(tree, hists, w, sysErrSquaredHists, err, errorOpt, additionalCuts, formulas, debugLevel)
        sumweight += w
        count += 1
    #propagate systematics to each histogram
    addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel)
    print "Sum of weights for this sample:",sumweight
    return sumweight

def loopTrees(treeDict, weightF, cuts="", hists={}, weightHists={}, sfHists={}, scale=1.0, weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], errorOpt=None, fillF=basicFill, sfVars=("MR","Rsq"), sysVars=("MR","Rsq"), boxName="NONE", auxSFs={}, debugLevel=0):
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
        auxSFHists = {name:sfHists[name] for name in auxSFs} #for misc reweightings
        sumweights += loopTree(treeDict[name], weightF, cuts, hists[name], weightHists, sfHistToUse, scale, fillF, sfVars, sysVars, weightOpts, errorOpt, process=name+"_"+boxName, auxSFs=auxSFs, auxSFHists=auxSFHists, debugLevel=debugLevel)
    print "Sum of event weights for all processes:",sumweights

