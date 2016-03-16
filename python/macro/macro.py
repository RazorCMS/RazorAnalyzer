import os
import ROOT as rt
import copy
import math
from array import array
from ast import literal_eval

#local imports
from PlotFit import setFFColors
from RunCombine import exec_me
from razorWeights import applyMTUncertainty1D, applyMTUncertainty2D, getSFHistNameForErrorOpt
from printJson import walk
import plotting

def walkDictionary(d, path=""):
    """
    Recursively iterate over a nested dictionary.
    Returns ordered pairs (p, x) where x is a value stored at some level of the dictionary and 
    p is the concatenation of all keys (joined with '/') needed to reach x from the top level.
    """
    for key,val in d.iteritems():
        #add key to path
        newpath = path+'/'+str(key)
        if newpath.startswith('/'): 
            newpath = newpath[1:]
        #yield or recurse
        if isinstance(val, dict):
            for x in walkDictionary(val, newpath):
                yield x
        else:
            yield (newpath, val)

def exportHists(hists, outFileName='hists.root', outDir='.', useDirectoryStructure=True, varName=None, delete=True, debugLevel=0):
    """
    Writes histograms from the given dictionary to a ROOT file.  
    If useDirectoryStructure is True, histograms will be put in folders according to the structure of the dictionary.
    If varName is provided, only histograms of the provided variable will be saved.
    """

    print "Storing histograms in file",outFileName
    outFile = rt.TFile(outDir+'/'+outFileName, 'recreate')
    for pair in walkDictionary(hists): #get (path_to_histogram, histogram) pairs
        if useDirectoryStructure:
            path = pair[0]
            if varName is not None:
                if str(varName) not in path: #skip histograms of other variables
                    continue
                else: #remove variable name from path
                    path = path.replace(str(varName),'').replace('//','/')
            tdir = outFile.GetDirectory(path)
            if tdir==None:
                print "Making directory",path
                outFile.mkdir(path)
                tdir = outFile.GetDirectory(path)
                tdir.cd()
                print "Writing histogram",pair[1].GetName(),"to directory",path
                pair[1].Write()
                if delete: pair[1].Delete()
        else:
            print "Writing histogram",pair[1].GetName()
            pair[1].Write()
            if delete: pair[1].Delete()
    outFile.Close()

def importHists(inFileName='hists.root', debugLevel=0):
    """
    Creates a dictionary according to the directory structure of the input file, 
    and populates the dictionary with the histograms in the file.
    """

    hists = {}
    print "\nGetting histograms from file",inFileName
    inFile = rt.TFile.Open(inFileName)
    for dirpath, dirnames, filenames, tdirectory in walk(inFile):
        if len(filenames) > 0: #there are objects to retrieve
            dirInFile = dirpath.split(':')[-1].split('/')[1:] #get path to histogram as a list of subdirectories
            currentLayer = hists #keep track of current layer of the dictionary
            storedHist = False
            if dirInFile != ['']: #dirInFile will be [''] if there are no subdirectories
                for i,subdir in enumerate(dirInFile): #descend along path to histogram
                    if subdir not in currentLayer: 
                        if len(filenames) == 1 and i == len(dirInFile)-1:
                            #if there is only one histogram, store it now
                            try:
                                lastSubdir = literal_eval(subdir)
                                if debugLevel > 0:
                                    print "Converted",subdir,"to tuple"
                            except ValueError:
                                lastSubdir = subdir
                            currentLayer[lastSubdir] = tdirectory.Get(filenames[0])
                            currentLayer[lastSubdir].SetDirectory(0)
                            if debugLevel > 0: 
                                print "Retrieved histogram",filenames[0]
                            storedHist = True
                            break
                        currentLayer[subdir] = {} #create dict layer if needed
                    currentLayer = currentLayer[subdir] #now we go one level lower and repeat
            #now currentLayer is the lowest-level dictionary, where the remaining histograms should be inserted
            if not storedHist:
                for tkey in filenames:
                    currentLayer[tkey] = tdirectory.Get(tkey)
                    currentLayer[tkey].SetDirectory(0)
                    if debugLevel > 0: 
                        print "Retrieved histogram",tkey
    inFile.Close()
    return hists

def stitch(th1s = []):
    """Join the input histograms end-to-end to form one long 1D histogram.
    Bins in the output histogram have width 1 and the x-axis starts at 0."""
    if len(th1s) == 0:
        return
    #get total number of bins across all inputs
    totalNBins = sum([th1.GetNbinsX() for th1 in th1s])
    #concatenate hist names
    name = 'and'.join([th1.GetName() for th1 in th1s])
    #make output hist
    out = rt.TH1F(name, name, totalNBins, 0, totalNBins)
    out.SetDirectory(0)
    bn = 1
    for i,th1 in enumerate(th1s):
        for bx in range(1, th1.GetNbinsX()+1):
            out.SetBinContent(bn, th1.GetBinContent(bx))
            out.SetBinError(bn, th1.GetBinError(bx))
            bn += 1
    return out

def getBinBoundariesFromColumns(xbins, cols):
    """
    Get lower and upper bin edges in the x and y directions
    (Useful for making yield tables)
    """
    xLow = []
    xHigh = []
    yLow = []
    yHigh = []
    for i in range(len(xbins)-1):
        for j in range(len(cols[i])-1):
            xLow.append(xbins[i])
            xHigh.append(xbins[i+1])
            yLow.append(cols[i][j])
            yHigh.append(cols[i][j+1])
    return xLow, xHigh, yLow, yHigh

def makeTH2PolyFromColumns(name, title, xbins, cols):
    """
    Makes a TH2Poly histogram with rectangular bins. 
    xbins: list of bin edges in the x-direction.
    cols: 2D list indicating bin edges in the y-direction for each column
    """
    #construct TH2Poly
    poly = rt.TH2Poly(name, title, xbins[0], xbins[-1], cols[0][0], cols[0][-1])
    #add bins in each column
    for i in range(len(xbins)-1):
        for j in range(len(cols[i])-1):
            poly.AddBin( xbins[i], cols[i][j], xbins[i+1], cols[i][j+1] )

    poly.SetDirectory(0)
    poly.Sumw2()
    return poly

def fillTH2PolyFromTH2(th2, poly):
    """Fills the given TH2Poly histogram using the bin contents of the given TH2.  If multiple bins in the TH2 map to the same bin in the TH2Poly, the bin contents will be added and the errors will be summed in quadrature."""
    for bx in range(1, th2.GetNbinsX()+1):
        for by in range(1, th2.GetNbinsY()+1):
            #get bin number in TH2Poly
            centerX = th2.GetXaxis().GetBinCenter(bx)
            centerY = th2.GetYaxis().GetBinCenter(by)
            bn = poly.FindBin(centerX, centerY)
            poly.SetBinContent(bn, poly.GetBinContent(bn) + th2.GetBinContent(bx,by))
            #NOTE: TH2Poly has a bug that causes SetBinError(x) to set the error of bin x+1. Beware!!!
            poly.SetBinError(bn-1, ( poly.GetBinError(bn)**2 + th2.GetBinError(bx,by)**2 )**(0.5) )

def makeTH2PolyRatioHist(num, denom, xbins, cols):
    """Makes a TH2Poly using makeTH2PolyFromColumns on the provided numerator and denominator histograms; 
    then divides the two histograms and propagates uncertainty assuming uncorrelated num. and denom."""

    #make TH2Poly with custom binning
    numPoly = makeTH2PolyFromColumns( num.GetName()+"Poly", num.GetTitle()+"Poly", xbins, cols)
    denomPoly = makeTH2PolyFromColumns( denom.GetName()+"Poly", denom.GetTitle()+"Poly", xbins, cols)

    #populate with numerator and denominator histogram bin contents
    fillTH2PolyFromTH2(num, numPoly) 
    fillTH2PolyFromTH2(denom, denomPoly)

    #need to divide by hand
    for bn in range(1, numPoly.GetNumberOfBins()+1):
        if denomPoly.GetBinContent(bn) != 0:
            numPoly.SetBinContent(bn, numPoly.GetBinContent(bn)/denomPoly.GetBinContent(bn))
            #NOTE: TH2Poly has a bug that causes SetBinError(x) to set the error of bin x+1. Beware!!!
            numPoly.SetBinError(bn-1, ( ( numPoly.GetBinError(bn)/denomPoly.GetBinContent(bn) )**2 + ( denomPoly.GetBinError(bn)*numPoly.GetBinContent(bn)/(denomPoly.GetBinContent(bn))**2 )**2 )**(0.5) )

    return numPoly

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
                        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                            if name not in shape[1]: continue
                            curShape = shape[0]
                        else:
                            curShape = shape
                        if curShape+"Down" not in shapeHists[name]: shapeHists[name][curShape+"Down"] = {}
                        if curShape+"Up" not in shapeHists[name]: shapeHists[name][curShape+"Up"] = {}
                        shapeHists[name][curShape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Down")
                        shapeHists[name][curShape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Up")
            elif len(var) == 2: 
                #2D histograms
                title = [titles[v] if v in titles else v for v in var]
                hists[name][var] = rt.TH2F(regionName+var[0]+var[1]+name, ';'+title[0]+';'+title[1], len(bins[var[0]])-1, array('d',bins[var[0]]), len(bins[var[1]])-1, array('d',bins[var[1]]))
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                            if name not in shape[1]: continue
                            curShape = shape[0]
                        else:
                            curShape = shape
                        if curShape+"Down" not in shapeHists[name]: shapeHists[name][curShape+"Down"] = {}
                        if curShape+"Up" not in shapeHists[name]: shapeHists[name][curShape+"Up"] = {}
                        shapeHists[name][curShape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Down")
                        shapeHists[name][curShape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Up")
            elif len(var) == 3:
                #3D histograms
                title = [titles[v] if v in titles else v for v in var]
                hists[name][var] = rt.TH3F(regionName+var[0]+var[1]+var[2]+name, ';'+title[0]+';'+title[1]+';'+title[2], len(bins[var[0]])-1, array('d',bins[var[0]]), len(bins[var[1]])-1, array('d',bins[var[1]]), len(bins[var[2]])-1, array('d',bins[var[2]]))
                if samples is not None and name in samples:
                    for shape in shapeErrors:
                        if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                            if name not in shape[1]: continue
                            curShape = shape[0]
                        else:
                            curShape = shape
                        if curShape+"Down" not in shapeHists[name]: shapeHists[name][curShape+"Down"] = {}
                        if curShape+"Up" not in shapeHists[name]: shapeHists[name][curShape+"Up"] = {}
                        shapeHists[name][curShape+"Down"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Down")
                        shapeHists[name][curShape+"Up"][var] = hists[name][var].Clone(hists[name][var].GetName()+curShape+"Up")
        for var in hists[name]: 
            hists[name][var].Sumw2()
            hists[name][var].SetDirectory(0)
            if samples is not None and name in samples:
                for shape in shapeErrors:
                    if not isinstance(shape,basestring):
                        if name not in shape[1]: continue
                        curShape = shape[0]
                    else:
                        curShape = shape
                    shapeHists[name][curShape+"Down"][var].Sumw2()
                    shapeHists[name][curShape+"Up"][var].Sumw2()
                    shapeHists[name][curShape+"Down"][var].SetDirectory(0)
                    shapeHists[name][curShape+"Up"][var].SetDirectory(0)
    #deal correctly with Poisson errors on data
    if dataName in hists:
        for var in hists[dataName]:
            hists[dataName][var].SetBinErrorOption(rt.TH1.kPoisson)
    return hists,shapeHists

def propagateShapeSystematics(hists, samples, varList, shapeHists, shapeErrors, miscErrors=[], boxName="", debugLevel=0):
    """For each bin of the central histogram, add the appropriate uncertainties in quadrature with the statistical uncertainty.
    List of arguments is similar to razorMacros.makeControlSampleHists
    """
    for var in varList:
        for name in samples:
            for shape in shapeErrors:
                if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                    if name not in shape[1]: continue
                    curShape = shape[0]
                else:
                    curShape = shape
                if debugLevel > 0: print "Adding",curShape,"uncertainty in quadrature with",name,"errors for",var
                #loop over histogram bins
                for bx in range(hists[name][var].GetSize()+1):
                    #use difference between Up and Down histograms as uncertainty
                    sysErr = abs(shapeHists[name][curShape+'Up'][var].GetBinContent(bx) - shapeHists[name][curShape+'Down'][var].GetBinContent(bx))/2.0
                    #add in quadrature with existing error
                    oldErr = hists[name][var].GetBinError(bx)
                    hists[name][var].SetBinError(bx, (oldErr**2 + sysErr**2)**(0.5))
                    if debugLevel > 0 and sysErr > 0: print curShape,": Error on bin ",bx,"increases from",oldErr,"to",hists[name][var].GetBinError(bx),"after adding",sysErr,"in quadrature"
            for source in miscErrors:
                #MT uncertainty (deprecated)
                if source.lower() == "mt" and var == "MR":
                    if isinstance(var, basestring): #1D histogram
                        applyMTUncertainty1D(hists[name][var], process=name+"_"+boxName, debugLevel=debugLevel)
                    else: #2D histogram
                        applyMTUncertainty2D(hists[name][var], process=name+"_"+boxName, debugLevel=debugLevel)

def subtractBkgsInData(process, hists={}, dataName="Data", debugLevel=0):
    """
    Subtracts all MC backgrounds (except 'process') from each data yield histogram.
    The data histogram is modified in-place in the hists dictionary
    """
    #warn if the needed histograms are not found
    if dataName not in hists:
        print "Error in subtractBkgsInData: target data histogram (",dataName,") was not found!"
        return
    bgProcesses = [mcProcess for mcProcess in hists if mcProcess != process and mcProcess != dataName and mcProcess != "Fit"] #get list of non-data, non-fit, non-signal samples
    if debugLevel > 0:
        print "Will isolate",process,"process by subtracting these backgrounds from the data yield:"
        print bgProcesses
    #subtract backgrounds in each data histogram
    for var in hists[dataName]:
        for p in bgProcesses:
            #make sure relevant background histogram exists
            if var not in hists[p]:
                print "Warning in subtractBkgsInData: could not find",var," in hists[",p,"]."
                continue
            #subtract it
            if debugLevel > 0: 
                print "Subtracting",p,"from",dataName,"distribution for",var
            hists[dataName][var].Add(hists[p][var], -1) 
        #zero any negative bins
        for bx in range(1, hists[dataName][var].GetSize()+1):
            if hists[dataName][var].GetBinContent(bx) < 0:
                print "Zeroing negative prediction",hists[dataName][var].GetBinContent(bx),"for",dataName
                hists[dataName][var].SetBinContent(bx, 0)

def basicPrint(histDict, mcNames, varList, c, printName="Hist", dataName="Data", logx=False, ymin=0.1, lumistr="40 pb^{-1}", boxName=None, btags=None, comment=True, blindBins=None, nsigmaFitData=None, nsigmaFitMC=None, doDensity=False, printdir=".", special="", unrollBins=(None,None), vartitles={}):
    """Make stacked plots of quantities of interest, with data overlaid"""
    #format MC histograms
    for name in mcNames: 
        for var in histDict[name]: plotting.setHistColor(histDict[name][var], name)
    titles = {name:name for name in mcNames}

    #get data histograms
    if dataName in histDict: dataHists = histDict[dataName]
    else: dataHists = None

    #make correct comment string
    if comment:
        if boxName is None or boxName == "":
            commentstr = printName
        else:
            commentstr = '#it{'+boxName
            if btags is not None and btags >= 0:
                commentstr += " "+str(btags)+" b-tag"
            #if 'sideband' in special:
                #commentstr += " Sideband Fit"
            #elif 'full' in special:
                #commentstr += ' Full Fit'
            commentstr += '}'
    else:
        commentstr = ""

    legend=None
    plotFit = ("Fit" in histDict)
    for i,var in enumerate(varList): 
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
            plotting.plot_basic_2D(c, mc=mcPrediction, data=obsData, fit=fitPrediction, xtitle=xtitle, ytitle=ytitle, printstr=var[0]+var[1]+printName, lumistr=lumistr, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, nsigmaFitData=nsigmaFitData, nsigmaFitMC=nsigmaFitMC, mcDict=mcDict, mcSamples=mcNames, ymin=ymin, unrollBins=unrollBins, printdir=printdir)
            #do MC total (no stack)
            plotting.plot_basic_2D(c, mc=mcPrediction, data=obsData, fit=fitPrediction, xtitle=xtitle, ytitle=ytitle, printstr=var[0]+var[1]+printName+'MCTotal', lumistr=lumistr, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, nsigmaFitData=nsigmaFitData, nsigmaFitMC=nsigmaFitMC, ymin=ymin, unrollBins=unrollBins, printdir=printdir)
            #draw bin mapping
            if unrollBins[0] is not None and unrollBins[1] is not None:
                plotting.drawUnrolledBinMapping(c, unrollBins, xtitle=xtitle, ytitle=ytitle, printstr=var[0]+var[1]+printName+"BINNING", printdir=printdir)
            #print prediction in each bin
            if obsData is not None and obsData != 0:
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
        stack = plotting.makeStack(varHists, mcNames, var)
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
            legend = plotting.makeLegend(varHists, titles, reversed(mcNames))
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
        if var in ['MR','Rsq','MR_NoW',"Rsq_NoW","MR_NoZ","Rsq_NoZ", "lep1.Pt()", "genlep1.Pt()","leadingGenLeptonPt"]: logx = True
        else: logx = False
        if blindBins is None:
            plotting.plot_basic(c, mc=stack, data=obsData, fit=fitPrediction, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=var+"_"+printName, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, printdir=printdir)
        else:
            plotting.plot_basic(c, mc=stack, data=None, fit=fitPrediction, leg=legend, xtitle=xtitle, ytitle=ytitle, printstr=var+"_"+printName, logx=logx, lumistr=lumistr, ymin=ymin, commentstr=commentstr, saveroot=True, savepdf=True, savepng=True, printdir=printdir)
        legend.Delete()

def transformVarsInString(string, varNames, suffix):
    outstring = copy.copy(string)
    for var in varNames:
        outstring = outstring.replace(var, var+suffix) 
    return outstring

JETVARS = ["MR","Rsq","nBTaggedJets","dPhiRazor","leadingJetPt","subleadingJetPt","nSelectedJets","nJets80","mT","box"] #quantities susceptible to jet uncertainties
def transformVarString(event, string, errorOpt, process="", debugLevel=0):
    outstring = string
    if errorOpt == "jesUp":
        outstring = transformVarsInString(string, JETVARS, "_JESUp")
    elif errorOpt == "jesDown":
        outstring = transformVarsInString(string, JETVARS, "_JESDown")
    elif errorOpt == "mesUp":
        outstring = transformVarsInString(string, JETVARS, "_MESUp")
    elif errorOpt == "mesDown":
        outstring = transformVarsInString(string, JETVARS, "_MESDown")
    elif errorOpt == "eesUp":
        outstring = transformVarsInString(string, JETVARS, "_EESUp")
    elif errorOpt == "eesDown":
        outstring = transformVarsInString(string, JETVARS, "_EESDown")
    elif errorOpt == "jerUp":
        outstring = transformVarsInString(string, JETVARS, "_JERUp")
    elif errorOpt == "jerDown":
        outstring = transformVarsInString(string, JETVARS, "_JERDown")

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
                sysErrSquaredHists[varName].Fill(varValue, sysErrSquared)
                if debugLevel > 1: print "Sys. Error =",sysErr,"; Filling (w*sysErr)^2 histogram with",sysErrSquared
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
                if debugLevel > 1: print "Sys. Error =",sysErr,"; Filling (w*sysErr)^2 histogram with",sysErrSquared

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
    if debugLevel > 1:
        print "Getting scale factor from histogram",sfHist.GetName()
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
    #TH2Poly case
    if sfHist.InheritsFrom('TH2Poly'):
        scaleFactor = sfHist.GetBinContent(sfHist.FindBin(*var))
        scaleFactorErr = sfHist.GetBinError(sfHist.FindBin(*var))
    #TH2F case
    else:
        scaleFactor = sfHist.GetBinContent(sfHist.FindFixBin(*var))
        scaleFactorErr = sfHist.GetBinError(sfHist.FindFixBin(*var))

    #add protection against unphysics scale factors
    #if scale factor is 0, then use scale factor of 1 and add 100% systematic uncertainty
    if scaleFactor == 0:
        scaleFactor = 1.0
        scaleFactorErr = math.sqrt( 1 + scaleFactorErr*scaleFactorErr )
    #if scale factor is > 2.0, then use scale factor of 1 and add difference between scale factor and 1.0 as a systematic uncertainty
    if scaleFactor > 2.0:
        scaleFactorErr = math.sqrt( (scaleFactor-1)*(scaleFactor-1) + scaleFactorErr*scaleFactorErr )
        scaleFactor = 1.0        
    #protect against NaN uncertainties
    if math.isnan(scaleFactorErr):
        scaleFactorErr = 0.0

    if debugLevel > 1: print "Applying scale factor: ",scaleFactor,"with error",scaleFactorErr
    return (scaleFactor, scaleFactorErr)

def addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel=0):
    """For each histogram in hists, look for the corresponding histogram in sysErrSquaredHists.
    Treats the values of sysErrSquaredHists as sums of (weight*error)^2 in each bin, and adds these errors in quadrature with the existing bin errors in hists"""
    for name in hists:
        if name in sysErrSquaredHists:
            if debugLevel > 0: print "Including systematic errors on ",name
            for bx in range(hists[name].GetSize()+1):
                squaredError = sysErrSquaredHists[name].GetBinContent(bx)
                oldErr = hists[name].GetBinError(bx)
                hists[name].SetBinError(bx,(oldErr*oldErr + squaredError)**(0.5))
                if debugLevel > 0 and squaredError > 0: print name,": Error on bin",bx,"increases from",oldErr,"to",hists[name].GetBinError(bx),"after adding",(squaredError**(0.5)),"in quadrature"

def propagateScaleFactorStatErrors(sysErrSquaredHists, upHists, downHists, debugLevel=0):
    """Use scale factor uncertainties to make up and down shape histograms"""
    for var in upHists:
        if var in sysErrSquaredHists:
            if debugLevel > 0: print "Making shape histogram for scale factor errors on",var
            for bx in range(1,sysErrSquaredHists[var].GetSize()+1):
                squaredError = sysErrSquaredHists[var].GetBinContent(bx)
                #increase up histogram bin content by 1 sigma, decrease down histogram bin content
                centralValue = upHists[var].GetBinContent(bx)
                upHists[var].SetBinContent(bx, centralValue + (squaredError)**(0.5))
                if centralValue > 0:
                    percentChange = (upHists[var].GetBinContent(bx) - centralValue)/centralValue
                    downHists[var].SetBinContent(bx, centralValue/(1+percentChange))

def loopTree(tree, weightF, cuts="", hists={}, weightHists={}, sfHist=None, scale=1.0, fillF=basicFill, sfVars=("MR","Rsq"), statErrOnly=False, weightOpts=[], errorOpt=None, process="", auxSFs={}, auxSFHists={}, shapeHists={}, shapeNames=[], shapeSFHists={}, shapeAuxSFs={}, shapeAuxSFHists={}, propagateScaleFactorErrs=True, noFill=False, debugLevel=0):
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
        #add the reweighting variable to the formula list
        if pair[0] not in formulas:
            formulas[pair[0]] = rt.TTreeFormula(pair[0], pair[0], tree)
            formulas[pair[0]].GetNdata()

    #make TTreeFormulas for shape histogram scale factors
    #A list of TTreeFormulas is maintained in shapeAuxSFForms, and shapeAuxSFFormsLookup matches each shape uncertainty with the index of a TTreeFormula in shapeAuxSFForms.
    shapeAuxSFForms = []
    shapeAuxSFFormsLookup = {} #stores the number of the TTreeFormula for each shape uncertainty cut
    cutStrings = []
    for n in shapeNames:
        shapeAuxSFFormsLookup[n+'Up'] = {}
        shapeAuxSFFormsLookup[n+'Down'] = {}
        for name,pair in shapeAuxSFs[n+'Up'].iteritems():
            #add reweighting variable to formula list
            if pair[0] not in formulas and isinstance(pair[0], basestring):
                formulas[pair[0]] = rt.TTreeFormula(pair[0], pair[0], tree)
                formulas[pair[0]].GetNdata()
            if pair[1] not in cutStrings:
                #make the TTreeFormula for this cut string
                if debugLevel > 0: print "Making TTreeFormula for",name,"with formula",pair[1],"for reweighting",pair[0],"(error option",n+"Up)"
                shapeAuxSFFormsLookup[n+'Up'][name] = len(shapeAuxSFForms) #index of TTreeFormula to use
                shapeAuxSFForms.append(rt.TTreeFormula("shapeAuxSFCuts"+str(len(shapeAuxSFForms)), pair[1], tree))
                shapeAuxSFForms[-1].GetNdata()
                cutStrings.append(pair[1])
            else:  
                #refer to the matching TTreeFormula that already exists
                index = cutStrings.index(pair[1])
                shapeAuxSFFormsLookup[n+'Up'][name] = index
        for name,pair in shapeAuxSFs[n+'Down'].iteritems():
            #add reweighting variable to formula list
            if pair[0] not in formulas and isinstance(pair[0], basestring):
                formulas[pair[0]] = rt.TTreeFormula(pair[0], pair[0], tree)
                formulas[pair[0]].GetNdata()
            if pair[1] not in cutStrings:
                #make the TTreeFormula for this cut string
                if debugLevel > 0: print "Making TTreeFormula for",name,"with formula",pair[1],"for reweighting",pair[0],"(error option",n+"Down)"
                shapeAuxSFFormsLookup[n+'Down'][name] = len(shapeAuxSFForms) #index of TTreeFormula to use
                shapeAuxSFForms.append(rt.TTreeFormula("shapeAuxSFCuts"+str(len(shapeAuxSFForms)), pair[1], tree))
                shapeAuxSFForms[-1].GetNdata()
                cutStrings.append(pair[1])
            else:  
                #refer to the matching TTreeFormula that already exists
                index = cutStrings.index(pair[1])
                shapeAuxSFFormsLookup[n+'Down'][name] = index

    #transform cuts 
    if errorOpt is not None:
        if debugLevel > 0: print "Error option is:",errorOpt
        cuts = transformVarString(tree, cuts, errorOpt, process=process, debugLevel=debugLevel+1)
    #create histograms for scale factor systematics
    sysErrSquaredHists = {}
    if not statErrOnly:
        for name in hists: 
            sysErrSquaredHists[name] = hists[name].Clone(hists[name].GetName()+"SFERRORS")
            sysErrSquaredHists[name].Reset()
            if debugLevel > 0: print "Created temp histogram",sysErrSquaredHists[name].GetName(),"to hold systematic errors from",sfVars,"scale factors"
    count = 0
    sumweight = 0.0
    if not noFill:
        #get list of entries passing the cuts
        tree.Draw('>>elist', cuts, 'entrylist')
        elist = rt.gDirectory.Get('elist')
        print "Total entries passing cuts:",elist.GetN()
        while True:
            #load the next entry
            entry = elist.Next()
            if entry == -1: break
            if count > 0 and count % 100000 == 0: print "Processing entry",count
            elif debugLevel > 0 and count % 10000 == 0: print "Processing entry",count
            elif debugLevel > 1: print "Processing entry",count
            tree.GetEntry(entry)
            w = weightF(tree, weightHists, scale, weightOpts, errorOpt, debugLevel=debugLevel)

            err = 0.0
            #get scale factor
            if sfHist is not None: 
                sf, err = getScaleFactorAndError(tree, sfHist, sfVars, formulas, debugLevel=debugLevel)
                #apply scale factor to event weight
                w *= sf
            for name in auxSFs: #apply misc scale factors (e.g. veto lepton correction)
                if auxSFForms[name].EvalInstance(): #check if this event should be reweighted
                    auxSF, auxErr = getScaleFactorAndError(tree, auxSFHists[name], sfVars=auxSFs[name][0], formulas=formulas, debugLevel=debugLevel)
                    w *= auxSF
                    err = (err*err + auxErr*auxErr)**(0.5)
            #protection for case of infinite-weight events
            if math.isinf(w):
                if debugLevel > 0: 
                    print "Warning: infinite-weight event encountered!"
                continue
            if math.isnan(err):
                print "Error: err is nan!"
            #fill with weight
            fillF(tree, hists, w, sysErrSquaredHists, err, errorOpt, additionalCuts=None, formulas=formulas, debugLevel=debugLevel)

            ###PROCESS EXTRA SHAPE HISTOGRAMS
            #get weights for filling additional shape histograms, if provided
            wsUp = [weightF(tree, weightHists, scale, weightOpts, e+'Up', debugLevel=debugLevel) for e in shapeNames]
            wsDown = [weightF(tree, weightHists, scale, weightOpts, e+'Down', debugLevel=debugLevel) for e in shapeNames]
            #evaluate all TTreeFormulas 
            shapeAuxSFFormResults = [form.EvalInstance() for form in shapeAuxSFForms]
            if debugLevel > 1 and len(shapeAuxSFFormResults) > 0:
                print "Cut decisions for systematic uncertainty scale factors:",shapeAuxSFFormResults
            #apply scale factor to shape up/down weights too (if needed)
            for i,n in enumerate(shapeNames):
                if debugLevel > 1: print "Error option:",n
                #main scale factor
                if shapeSFHists[n+'Up'] is not None:
                    shapeSFUp, shapeErrUp = getScaleFactorAndError(tree, shapeSFHists[n+'Up'], sfVars, formulas, debugLevel=debugLevel)
                    wsUp[i] *= shapeSFUp
                if shapeSFHists[n+'Down'] is not None:
                    shapeSFDown, shapeErrDown = getScaleFactorAndError(tree, shapeSFHists[n+'Down'], sfVars, formulas, debugLevel=debugLevel)
                    wsDown[i] *= shapeSFDown
                #misc scale factors
                for name in shapeAuxSFs[n+'Up']:
                    if shapeAuxSFFormResults[ shapeAuxSFFormsLookup[n+'Up'][name] ]:
                        shapeAuxSFUp, shapeAuxErrUp = getScaleFactorAndError(tree, shapeAuxSFHists[n+'Up'][name], sfVars=shapeAuxSFs[n+'Up'][name][0], formulas=formulas, debugLevel=debugLevel)
                        wsUp[i] *= shapeAuxSFUp
                for name in shapeAuxSFs[n+'Down']:
                    if shapeAuxSFFormResults[ shapeAuxSFFormsLookup[n+'Down'][name] ]:
                        shapeAuxSFDown, shapeAuxErrDown = getScaleFactorAndError(tree, shapeAuxSFHists[n+'Down'][name], sfVars=shapeAuxSFs[n+'Down'][name][0], formulas=formulas, debugLevel=debugLevel)
                        wsDown[i] *= shapeAuxSFDown
                #fill up/down shape histograms
                fillF(tree, shapeHists[n+'Up'], wsUp[i], formulas=formulas, debugLevel=debugLevel)
                fillF(tree, shapeHists[n+'Down'], wsDown[i], formulas=formulas, debugLevel=debugLevel)
            #####

            sumweight += w
            count += 1
    #propagate systematics to each histogram
    if 'sfstatttjets' in shapeNames:
        propagateScaleFactorStatErrors(sysErrSquaredHists, upHists=shapeHists['sfstatttjetsUp'], downHists=shapeHists['sfstatttjetsDown'], debugLevel=debugLevel)
    elif 'sfstatwjets' in shapeNames:
        propagateScaleFactorStatErrors(sysErrSquaredHists, upHists=shapeHists['sfstatwjetsUp'], downHists=shapeHists['sfstatwjetsDown'], debugLevel=debugLevel)
    elif 'sfstatzinv' in shapeNames:
        propagateScaleFactorStatErrors(sysErrSquaredHists, upHists=shapeHists['sfstatzinvUp'], downHists=shapeHists['sfstatzinvDown'], debugLevel=debugLevel)
    elif propagateScaleFactorErrs:
        addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel)
    print "Sum of weights for this sample:",sumweight
    return sumweight

def loopTrees(treeDict, weightF, cuts="", hists={}, weightHists={}, sfHists={}, scale=1.0, weightOpts=[], errorOpt=None, fillF=basicFill, sfVars=("MR","Rsq"), statErrOnly=False, boxName="NONE", auxSFs={}, shapeHists={}, shapeNames=[], shapeAuxSFs={}, noFill=False, propagateScaleFactorErrs=True, debugLevel=0):
    """calls loopTree on each tree in the dictionary.  
    Here hists should be a dict of dicts, with hists[name] the collection of histograms to fill using treeDict[name]"""
    sumweights=0.0
    for name in treeDict: 
        if name not in hists: continue
        print("Filling histograms for tree "+name)

        #get correct scale factor histogram
        sfHistToUse = None
        if name in sfHists:
            sfHistToUse = sfHists[name]
            print("Using scale factors from histogram "+sfHistToUse.GetName())
        #get appropriate scale factor histograms for misc reweightings
        auxSFHists = {name:sfHists[name] for name in auxSFs} 
        #get correct variables for scale factor reweighting.
        #if sfVars is a dictionary, get the appropriate value from it.  otherwise use sfVars directly.
        sfVarsToUse = sfVars
        if sfHistToUse is not None and isinstance(sfVars,dict):
            sfVarsToUse = sfVars[name]
            print "Reweighting in",sfVarsToUse


        #prepare additional shape histograms if needed
        shapeHistsToUse = {}
        shapeNamesToUse = []
        shapeAuxSFsToUse = {}
        shapeAuxSFHists = {}
        shapeSFHists = {}
        for shape in shapeNames:
            if not isinstance(shape,basestring): #tuple (shape, [list of processes])
                if name not in shape[1]: continue
                curShape = shape[0]
            else:
                curShape = shape
            if curShape+'Up' in shapeHists[name]:
                shapeHistsToUse[curShape+'Up'] = shapeHists[name][curShape+'Up']
                shapeHistsToUse[curShape+'Down'] = shapeHists[name][curShape+'Down']
                shapeNamesToUse.append(curShape)
                shapeAuxSFsToUse[curShape+'Up'] = shapeAuxSFs[curShape+'Up']
                shapeAuxSFsToUse[curShape+'Down'] = shapeAuxSFs[curShape+'Down']
                shapeAuxSFHists[curShape+'Up'] = {n:sfHists[n] for n in shapeAuxSFsToUse[curShape+'Up']} 
                shapeAuxSFHists[curShape+'Down'] = {n:sfHists[n] for n in shapeAuxSFsToUse[curShape+'Down']} 
                shapeSFHists[curShape+'Up'] = None
                shapeSFHists[curShape+'Down'] = None
                if name in sfHists:
                    shapeSFHists[curShape+'Up'] = sfHists[getSFHistNameForErrorOpt(curShape+'Up', name)]
                    shapeSFHists[curShape+'Down'] = sfHists[getSFHistNameForErrorOpt(curShape+'Down', name)]
        if debugLevel > 0:
            print "Will fill histograms for these shape uncertainties:"
            print shapeNamesToUse

        sumweights += loopTree(treeDict[name], weightF, cuts, hists[name], weightHists, sfHistToUse, scale, fillF, sfVarsToUse, statErrOnly, weightOpts, errorOpt, process=name+"_"+boxName, auxSFs=auxSFs, auxSFHists=auxSFHists, shapeHists=shapeHistsToUse, shapeNames=shapeNamesToUse, shapeSFHists=shapeSFHists, shapeAuxSFs=shapeAuxSFsToUse, shapeAuxSFHists=shapeAuxSFHists, noFill=noFill, propagateScaleFactorErrs=propagateScaleFactorErrs, debugLevel=debugLevel)
    print "Sum of event weights for all processes:",sumweights

def correctScaleFactorUncertaintyForSignalContamination(centralHist, upHist, downHist, sfHist, contamHist, debugLevel=0):
    """
    sigHist should be the histogram that needs to be corrected.
    sfHist should be the histogram of scale factors.
    contamHist should be a histogram of the same binning as sfHist giving % signal contamination in each bin of sfHist.
    Increases the uncertainty on each bin of sfHist to account for the level of signal contamination.
    Assumes that the signal histogram uncertainties reflect MC statistics only!  
    (i.e. no other systematics have been propagated)
    Supports TH2s as well as TH2Polys for the scale factor histogram
    """
    if debugLevel > 0:
        print "Adding uncertainty for signal contamination in",sfHist.GetName()
    if sfHist.InheritsFrom('TH2Poly'):
        for bx in range(1,centralHist.GetNbinsX()+1):
            for by in range(1, centralHist.GetNbinsY()+1):
                #find the scale factor bin corresponding to this signal region bin
                xCoord = centralHist.GetXaxis().GetBinCenter(bx)
                yCoord = centralHist.GetYaxis().GetBinCenter(by)
                sfBin = sfHist.FindBin(xCoord, yCoord)
            
                #get error, scale factor error, and level of signal contamination
                sf = sfHist.GetBinContent(sfBin)
                contam = contamHist.GetBinContent(sfBin)
                contamErr = sf*contam

                curErr = upHist.GetBinContent(bx,by) - centralHist.GetBinContent(bx,by)
                newErr = ( curErr**2 + contamErr**2 )**(0.5)

                upHist.SetBinContent(bx,by, centralHist.GetBinContent(bx,by) + newErr)
                if centralHist.GetBinContent(bx,by) > 0:
                    percentChange = newErr/centralHist.GetBinContent(bx,by)
                    downHist.SetBinContent(bx,by, centralHist.GetBinContent(bx,by)/(1+percentChange))

                if debugLevel > 0:
                    print "Signal contamination in bin",bx,by,"is",contam,"; uncertainty goes from",curErr,"to",newErr
    else:
        print "Error in correctScaleFactorUncertaintyForSignalContamination: function implemented for TH2Poly only!"
        sys.exit()
            
def getExcludedSignalStrength(dirName, model, mGluino=-1, mStop=-1, mLSP=-1, debugLevel=0): 
    """ Retrieve the expected signal exclusion for the given model, using previously computed limits """ 
    if mLSP < 0:
        print "Error in getExcludedSignalStrength: please specify mLSP!"
        return 0

    #open file and get limit results
    fName = model+'_MultiJet_EleMultiJet_MuMultiJet_results.root'
    hName = "xsecUL_Exp_%s_MultiJet_EleMultiJet_MuMultiJet"%(model)
    resultFile = rt.TFile.Open(dirName+'/'+fName)
    xsecULHist = resultFile.Get(hName)
    if not xsecULHist:
        print "Error in getExcludedSignalStrength: histogram",hName,"not found in",dirName+'/'+fName
        return 0

    #get excluded cross section
    if 'T2' in model:
        xsecUL = xsecULHist.GetBinContent(xsecULHist.FindFixBin(mStop,mLSP))
    else:
        xsecUL = xsecULHist.GetBinContent(xsecULHist.FindFixBin(mGluino,mLSP))

    #get theoretical cross section
    if mGluino!=-1:
        for line in open('data/gluino13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mGluino))==line.split(',')[0]:
                xsecTheory = float(line.split(',')[1]) #pb
    elif mStop!=-1:
        for line in open('data/stop13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mStop))==line.split(',')[0]:
                xsecTheory = float(line.split(',')[1]) #pb
    else:
        print "Error in getExcludedSignalStrength: please specify either mStop or mGluino!"
        return 0
    if debugLevel > 0:
        print "Expected UL",xsecUL,"pb"
        print "Theory cross section",xsecTheory,"pb"

    #return the ratio
    return xsecUL/xsecTheory
