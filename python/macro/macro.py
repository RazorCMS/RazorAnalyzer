import ROOT as rt
import copy
import array

def setFFColors(hNS, minZ=-5.1, maxZ=5.1):
    Red = array.array('d',  [0.00, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00])
    Green = array.array('d',[0.00, 0.70, 0.90, 1.00, 0.90, 0.70, 0.00])
    Blue = array.array('d', [1.00, 1.00, 1.00, 1.00, 0.90, 0.70, 0.00])
    Length =array.array('d',[0.00, 0.20, 0.35, 0.50, 0.65, 0.8, 1.00]) # colors get darker faster at 4sigma
    rt.TColor.CreateGradientColorTable(7,Length,Red,Green,Blue,999)
    hNS.SetMaximum(maxZ)
    hNS.SetMinimum(minZ) # so the binning is 0 2 4
    hNS.SetContour(999)

def basicPrint(histDict, mcNames, varList, c, printName="Hist", dataName="Data", logx=False, ymin=0.1, lumistr="40 pb^{-1}"):
    """Make stacked plots of quantities of interest, with data overlaid"""
    #format MC histograms
    for name in mcNames: 
        for var in histDict[name]: setHistColor(histDict[name][var], name)
    titles = {name:name for name in mcNames}

    #get data histograms
    dataHists = histDict[dataName]

    legend=None
    for i,var in enumerate(varList): 
        #for MR and Rsq, make 2D plots
        if var == ('MR','Rsq'):
            mcPrediction = histDict[mcNames[0]][var].Clone("mcPrediction")
            mcPrediction.Reset()
            for name in mcNames: 
                mcPrediction.Add(histDict[name][var])
            commentstr = printName+" Box"
            plot_basic_2D(c, mc=mcPrediction, data=dataHists[var], xtitle='MR', ytitle='Rsq', printstr='Razor_'+printName, lumistr=lumistr, commentstr=commentstr, saveroot=True)
        #for other variables make 1D plots
        if not isinstance(var, basestring): continue #only consider strings
        varHists = {name:histDict[name][var] for name in mcNames}
        if not legend:
            legend = makeLegend(varHists, titles, reversed(mcNames))
            legend.AddEntry(dataHists[var], dataName)
        stack = makeStack(varHists, mcNames, var)
        plot_basic(c, mc=stack, data=dataHists[var], leg=legend, xtitle=var, printstr=var+"_"+printName, logx=logx, lumistr=lumistr, ymin=ymin, saveroot=True)

def transformVarString(string, errorOpt, debugLevel=0):
    outstring = copy.copy(string)
    if errorOpt == "jesUp":
        outstring = outstring.replace("MR", "MR_JESUp")
        outstring = outstring.replace("Rsq", "Rsq_JESUp")
        outstring = outstring.replace("nBTaggedJets", "nBTaggedJets_JESUp")
        outstring = outstring.replace("dPhiRazor", "dPhiRazor_JESUp")
        outstring = outstring.replace("leadingJetPt", "leadingJetPt_JESUp")
        outstring = outstring.replace("subleadingJetPt", "subleadingJetPt_JESUp")
        outstring = outstring.replace("nSelectedJets", "nSelectedJets_JESUp")
        outstring = outstring.replace("nJets80", "nJets80_JESUp")
        outstring = outstring.replace("box", "box_JESUp")
    elif errorOpt == "jesDown":
        outstring = outstring.replace("MR", "MR_JESDown")
        outstring = outstring.replace("Rsq", "Rsq_JESDown")
        outstring = outstring.replace("nBTaggedJets", "nBTaggedJets_JESDown")
        outstring = outstring.replace("dPhiRazor", "dPhiRazor_JESDown")
        outstring = outstring.replace("leadingJetPt", "leadingJetPt_JESDown")
        outstring = outstring.replace("subleadingJetPt", "subleadingJetPt_JESDown")
        outstring = outstring.replace("nSelectedJets", "nSelectedJets_JESDown")
        outstring = outstring.replace("nJets80", "nJets80_JESDown")
        outstring = outstring.replace("box", "box_JESDown")
    if errorOpt == "jerUp":
        outstring = outstring.replace("MR", "MR_JERUp")
        outstring = outstring.replace("Rsq", "Rsq_JERUp")
        outstring = outstring.replace("nBTaggedJets", "nBTaggedJets_JERUp")
        outstring = outstring.replace("dPhiRazor", "dPhiRazor_JERUp")
        outstring = outstring.replace("leadingJetPt", "leadingJetPt_JERUp")
        outstring = outstring.replace("subleadingJetPt", "subleadingJetPt_JERUp")
        outstring = outstring.replace("nSelectedJets", "nSelectedJets_JERUp")
        outstring = outstring.replace("nJets80", "nJets80_JERUp")
        outstring = outstring.replace("box", "box_JERUp")
    elif errorOpt == "jerDown":
        outstring = outstring.replace("MR", "MR_JERDown")
        outstring = outstring.replace("Rsq", "Rsq_JERDown")
        outstring = outstring.replace("nBTaggedJets", "nBTaggedJets_JERDown")
        outstring = outstring.replace("dPhiRazor", "dPhiRazor_JERDown")
        outstring = outstring.replace("leadingJetPt", "leadingJetPt_JERDown")
        outstring = outstring.replace("subleadingJetPt", "subleadingJetPt_JERDown")
        outstring = outstring.replace("nSelectedJets", "nSelectedJets_JERDown")
        outstring = outstring.replace("nJets80", "nJets80_JERDown")
        outstring = outstring.replace("box", "box_JERDown")

    if debugLevel > 1:
        if outstring != string: print "For option",errorOpt,"Replacing string '",string,"' with '",outstring,"'"
    return outstring

def basicFill(tree, hists={}, weight=1.0, sysErrSquaredHists={}, sysErr=0.0, errorOpt=None, debugLevel=0):
    """Fills each histogram with the corresponding variable in the tree.
    'hists' should be a dictionary of histograms, with keys being the variable names to fill.
    Ex: hists['MR'] should be the histogram you want to fill with MR values.
    A key that is a tuple of variables (ex: ('MR','Rsq')) should be paired with a multidimensional histogram.
    In this case, the given variables will be used to fill the histogram."""
    for varName, hist in hists.iteritems(): 
        if isinstance(varName, basestring): #if varName is a string
            #transform variable name
            if errorOpt is not None: varName = transformVarString(varName, errorOpt, debugLevel=debugLevel)
            if debugLevel > 1: print "Filling",varName,"=",getattr(tree,varName),"with weight",weight
            hist.Fill(getattr(tree, varName), weight)
            if varName in sysErrSquaredHists: #for propagating systematic errors on the variables
                sysErrSquared = weight*weight*sysErr*sysErr
                sysErrSquaredHist[varName].Fill(getattr(tree, varName), sysErrSquared)
        else: #treat it as a tuple of variables that should be filled
            #transform each variable
            if errorOpt is not None: varName = tuple([transformVarString(v, errorOpt, debugLevel) for v in varName])
            toFill = [getattr(tree, v) for v in varName]+[weight]
            if debugLevel > 1: print "Filling",varName,":",toFill
            hist.Fill(*toFill)
            if varName in sysErrSquaredHists:
                sysErrSquared = weight*weight*sysErr*sysErr
                toFillErr = [getattr(tree, v) for v in varName]+[sysErrSquared]
                sysErrSquaredHists[varName].Fill(*toFillErr)

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
                    hists[name].SetBinError(bx,by,(hists[name].GetBinError(bx,by)**2 + squaredError)**(0.5))

def loopTree(tree, weightF, cuts="", hists={}, weightHists={}, sfHist=None, scale=1.0, fillF=basicFill, sfVars=("MR","Rsq"), sysVars=("MR", "Rsq"), weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], errorOpt=None, debugLevel=0):
    """Loop over a single tree and fill histograms.
    Returns the sum of the weights of selected events."""
    if debugLevel > 0: print ("Looping tree "+tree.GetName())
    if debugLevel > 0: print ("Cuts: "+cuts)
    #transform cuts 
    if errorOpt is not None:
        if debugLevel > 0: print "Error option is:",errorOpt
        cuts = transformVarString(cuts, errorOpt, debugLevel=debugLevel+1)
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
        err = 0.0
        if sfHist is not None: 
            sf, err = getScaleFactorAndError(tree, sfHist, sfVars, debugLevel)
            w *= sf
        fillF(tree, hists, w, sysErrSquaredHists, err, errorOpt, debugLevel)
        sumweight += w
        count += 1
    #propagate systematics to each histogram
    addToTH2ErrorsInQuadrature(hists, sysErrSquaredHists, debugLevel)
    print "Sum of weights for this sample:",sumweight
    return sumweight

def loopTrees(treeDict, weightF, cuts="", hists={}, weightHists={}, sfHists={}, scale=1.0, weightOpts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], errorOpt=None, fillF=basicFill, sfVars=("MR","Rsq"), sysVars=("MR","Rsq"), debugLevel=0):
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
        sumweights += loopTree(treeDict[name], weightF, cuts, hists[name], weightHists, sfHistToUse, scale, fillF, sfVars, sysVars, weightOpts, errorOpt, debugLevel)
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

colors = {"WJets":rt.kRed+1, "DYJets":rt.kBlue+1, "TTJets":rt.kGreen+2, "ZInv":rt.kCyan+1, "QCD":rt.kOrange+3, "SingleTop":rt.kOrange-3, "VV":rt.kViolet+3, "TTV":rt.kGreen-7}
def setHistColor(hist, name):
    """Sets histogram color"""
    if name in colors: hist.SetFillColor(colors[name])
    else: print "Warning in macro.py: histogram fill color not set for",name

def plot_basic(c, mc=0, data=0, fit=0, leg=0, xtitle="", ytitle="Number of events", ymin=None, ymax=None, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", ratiomin=0.5, ratiomax=1.5, saveroot=False, savepdf=False, savepng=True):
    """Plotting macro with options for data, MC, and fit histograms.  Creates data/MC ratio if able."""
    #setup
    c.Clear()
    c.cd()
    if data and mc: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.4, 1, 1)
    else: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.1, 1, 1)
    pad1.SetBottomMargin(0)
    pad1.SetLogx(logx)
    pad1.SetLogy(logy)
    pad1.Draw()
    pad1.cd()
    #draw MC
    if mc:
        mc.SetTitle("")
        mc.Draw("hist")
        if logy: mc.GetXaxis().SetMoreLogLabels()
        if not data: mc.GetXaxis().SetTitle(xtitle)
        mc.GetYaxis().SetTitle(ytitle)
        mc.GetYaxis().SetLabelSize(0.03)
        if data: mc.GetYaxis().SetTitleOffset(0.45)
        else: mc.GetYaxis().SetTitleOffset(0.50)
        mc.GetYaxis().SetTitleSize(0.05)
        if ymin is not None: mc.SetMinimum(ymin)
        if ymax is not None: mc.SetMaximum(ymax)
    #draw data
    if data:
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1)
        data.GetYaxis().SetTitle(ytitle)
        data.Draw("pesame")
    if fit:
        fitHist.SetLineWidth(2)
        fitHist.Draw("lsame")
    pad1.Modified()
    rt.gPad.Update()
    #make ratio data/MC
    if data and mc:
        histList = mc.GetHists()
        mcTotal = histList.First().Clone()
        mcTotal.Reset()
        for h in histList:
            mcTotal.Add(h)
        dataOverMC = data.Clone()
        dataOverMC.Divide(mcTotal)
        dataOverMC.SetTitle("")
        dataOverMC.GetXaxis().SetTitle(xtitle)
        dataOverMC.GetYaxis().SetTitle("Data / MC")
        dataOverMC.SetMinimum(ratiomin)
        dataOverMC.SetMaximum(ratiomax)
        dataOverMC.GetXaxis().SetLabelSize(0.1)
        dataOverMC.GetYaxis().SetLabelSize(0.08)
        dataOverMC.GetYaxis().SetTitleOffset(0.35)
        dataOverMC.GetXaxis().SetTitleOffset(1.00)
        dataOverMC.GetYaxis().SetTitleSize(0.08)
        dataOverMC.GetXaxis().SetTitleSize(0.08)
        dataOverMC.SetStats(0)
        if logx: dataOverMC.GetXaxis().SetMoreLogLabels()
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
    #draw data/MC
    c.cd()
    if data and mc: 
        pad2 = rt.TPad(printstr+"pad2",printstr+"pad2",0,0.0,1,0.4)
        pad2.SetTopMargin(0)
        pad2.SetTopMargin(0.008)
        pad2.SetBottomMargin(0.25)
        pad2.SetGridy()
        pad2.SetLogx(logx)
        pad2.Draw()
        pad2.cd()
        dataOverMC.Draw("pe")
        pad2.Modified()
        rt.gPad.Update()
    #save
    if savepng: c.Print(printstr+".png")
    if savepdf: c.Print(printstr+".pdf")
    if saveroot: c.Print(printstr+".root")

def draw2DHist(c, hist, xtitle="", ytitle="", ztitle="", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, drawErrs=False, palette=53, saveroot=False, savepdf=False, savepng=True):
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
    hist.GetZaxis().SetTitle(ztitle)
    hist.GetYaxis().SetLabelSize(0.03)
    hist.GetYaxis().SetTitleOffset(0.50)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.SetStats(0)
    if zmin is not None: hist.SetMinimum(zmin)
    if zmax is not None: hist.SetMaximum(zmax)
    hist.Draw("colz")
    if dotext:
        rt.gStyle.SetPaintTextFormat('4.1f')
        hist.SetMarkerSize(2.0)
        if not drawErrs: hist.Draw('textsame')
        else: hist.Draw('textesame')
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
        t3 = rt.TLatex(0.40, 0.84, commentstr)
        t3.SetNDC()
        t3.SetTextSize(0.04)
        t3.Draw()
    #save
    if savepng: c.Print(printstr+".png")
    if savepdf: c.Print(printstr+".pdf")
    if saveroot: c.Print(printstr+".root")

def make2DPullHistogram(h1, h2):
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
    return ret

def plot_basic_2D(c, mc=0, data=0, fit=0, xtitle="", ytitle="", ztitle="Number of events", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, saveroot=False, savepdf=False, savepng=True):
    """Plotting macro for data, MC, and/or fit yields.  Creates french flag plots comparing data/MC/fit if able."""
    #draw each histogram separately
    if mc:
        draw2DHist(c, mc, xtitle, ytitle, ztitle, zmin, zmax, printstr+'MC', lumistr=lumistr, commentstr=commentstr+", MC prediction")
        if data: #do (data - mc)/unc
            mcPulls = make2DPullHistogram(data,mc)
            draw2DHist(c, mcPulls, xtitle, ytitle, ztitle, None, None, printstr+'MCPulls', lumistr=lumistr, commentstr=commentstr+", (Data - MC)/#sigma", palette="FF", logz=False)
        if fit: #do (mc - fit)/unc
            mcFitPulls = make2DPullHistogram(mc,fit)
            draw2DHist(c, mcFitPulls, xtitle, ytitle, ztitle, None, None, printstr+'MCFitPulls', lumistr=lumistr, commentstr=commentstr+", (MC - Fit)/#sigma", palette="FF", logz=False)
    if data:
        draw2DHist(c, data, xtitle, ytitle, ztitle, zmin=max(0.1,zmin), printstr=printstr+'Data', lumistr=lumistr, commentstr=commentstr+", Data")
        if fit: #do (data - fit)/unc
            dataFitPulls = make2DPullHistogram(data,fit)
            draw2DHist(c, dataFitPulls, xtitle, ytitle, ztitle, None, None, printstr+'DataFitPulls', lumistr=lumistr, commentstr=commentstr+", (Data - Fit)/#sigma", palette="FF", logz=False)
    if fit:
        draw2DHist(c, fit, xtitle, ytitle, ztitle, zmin, zmax, printstr+'Fit', lumistr=lumistr, commentstr=commentstr+", Fit prediction")

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
        f.write('\\newgeometry{margin=0.2cm}\n')
        if landscape: f.write('\\begin{landscape}\n')
        f.write('\\begin{center}\n\\footnotesize\n\\begin{longtable}{|'+('|'.join(['c' for c in cols]))+'|}\n')
        f.write('\\caption{'+caption+'}\n\\endhead\n\\hline\n')
        f.write(' & '.join(headers)+' \\\\\n\\hline\n')
        for row in range(len(cols[0])):
            f.write((' & '.join([col[row] for col in cols]))+' \\\\\n\\hline\n')
        f.write('\\end{longtable}\n\\end{center}\n')
        if landscape: f.write('\\end{landscape}\n')
        f.write('\\restoregeometry\n')
        print "Created LaTeX scale factor table",(printstr+".tex")
