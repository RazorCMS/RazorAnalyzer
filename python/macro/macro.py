import ROOT as rt

def makeTreeDict(fileDict, treeName, debug=False):
    """Opens each file in fileDict, gets a tree called treeName, and returns a dict of trees"""
    trees = {}
    for name in fileDict:
        if debug: print("Loading tree "+treeName)
        trees[name] = fileDict[name].Get(treeName)
        if debug: print("Got tree containing "+str(trees[name].GetEntries())+" entries")
        assert trees[name]
    if debug: 
        print("Trees loaded: ") 
        print trees
    return trees

def loadWeightHists(debug=False):
    """Returns a dict with necessary reweighting histograms"""
    puWeightFileName = "TODO"
    wHists = {}
    #pileup weight histogram
    #(comment out until PU weights are available)
    #if debug: print("Opening pileup weight file "+puWeightFileName)
    #puFile = rt.TFile(puWeightFileName)
    #wHists["pileup"] = puFile.Get("TODO")
    #wHists["pileup"].SetDirectory(0)
    #assert wHists["pileup"]
    return wHists

def fillHist(hists, name, value, weight, debug):
    """Fill hists[name] with (value, weight)"""
    if name in hists:
        if debug: print("Filling histogram: "+name)
        hists[name].Fill(value, weight)
    else: 
        if debug: print("Warning in macro.py: histogram "+name+" not found!")

def loopTree(tree, cutF, weightF, weightHists, fillF, hists, maxEvents=-1, debug=False):
    """Loop over a single tree and fill histograms"""
    if debug: print ("Looping tree "+tree.GetName())
    for i,event in enumerate(tree):
        if maxEvents >= 0 and i >= maxEvents: break
        if i % 100000 == 0: print("Processing entry "+str(i))
        if not cutF(event, debug): continue
        w = weightF(event, weightHists, debug)
        fillF(event, hists, w, debug)

def loopTrees(treeDict, cutF, weightF, weightHists, fillF, histDict, maxEvents=-1, debug=False):
    """calls loopTree on each tree in the dictionary"""
    for name in treeDict: 
        loopTree(treeDict[name], cutF, weightF, weightHists, fillF, histDict[name], maxEvents, debug)

def makeStack(hists, ordering, title="Stack"):
    """Takes a dict of histograms and an ordered list of names, and returns a THStack containing the histograms stacked in the desired order"""
    stack = rt.THStack("thstack"+title.replace(" ",""), title)
    for name in ordering: 
        if name in hists:
            stack.Add(hists[name])
        else: 
            print("Warning in makeStack: histogram "+name+" not found in histogram dictionary!")
    return stack

def makeLegend(hists, names, ordering, x1=0.6, y1=0.6, x2=0.9, y2=0.9):
    """Takes a dict of histograms, a dict of histogram titles, and an ordered list of names, and returns a legend with the histograms in the desired order"""
    leg = rt.TLegend(x1, y1, x2, y2)
    for name in ordering: 
        leg.AddEntry(hists[name], names[name])
    return leg

def setHistColor(hist, name):
    """Sets histogram color according to the colors listed here"""
    colors = {}
    colors["WJets"] = 900
    colors["DYJets"] = 901
    colors["TTJets"] = 902
    colors["ZJetsNuNu"] = 903
    colors["QCD"] = 904
    colors["SingleTop"] = 905
    colors["VV"] = 906
    colors["TTV"] = 907
    red = rt.gROOT.GetColor(900)
    red.SetRGB(.42, .125, .125)
    blue = rt.gROOT.GetColor(901)
    blue.SetRGB(.106, .153, .282)
    green = rt.gROOT.GetColor(902)
    green.SetRGB(.102, .333, .102)
    lblue = rt.gROOT.GetColor(903)
    lblue.SetRGB(.267, .314, .455)
    orange = rt.gROOT.GetColor(904)
    orange.SetRGB(.42, .325, .125)
    lgreen = rt.gROOT.GetColor(905)
    lgreen.SetRGB(.29, .541, .29)
    dblue = rt.gROOT.GetColor(906)
    dblue.SetRGB(.012, .039, .114)
    dgreen = rt.gROOT.GetColor(907)
    dgreen.SetRGB(0., .137, 0.)

    if name in colors: hist.SetFillColor(colors[name])
    else: print("Warning in macro.py: histogram fill color not set")

def plot_basic(mc=0, data=0, fit=0, leg=0, xtitle="", ytitle="Number of events", ymin=0.1, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", ratiomin=0.5, ratiomax=1.5, saveroot=False, savepdf=False, savepng=True):
    """Plotting macro with options for data, MC, and fit histograms.  Creates data/MC ratio if able."""
    #setup
    c = rt.TCanvas("c", "c", 800, 600)
    c.Clear()
    c.cd()
    if data and mc: pad1 = rt.TPad("pad1", "pad1", 0, 0.4, 1, 1)
    else: pad1 = rt.TPad("pad1", "pad1", 0, 0.1, 1, 1)
    pad1.SetBottomMargin(0)
    pad1.SetLogx(logx)
    pad1.SetLogy(logy)
    pad1.Draw()
    pad1.cd()
    #draw MC
    if mc:
        mc.SetTitle("")
        mc.Draw("hist")
        if not data: mc.GetXaxis().SetTitle(xtitle)
        mc.GetYaxis().SetTitle(ytitle)
        mc.GetYaxis().SetLabelSize(0.03)
        if data: mc.GetYaxis().SetTitleOffset(0.45)
        else: mc.GetYaxis().SetTitleOffset(0.50)
        mc.GetYaxis().SetTitleSize(0.05)
        mc.SetMinimum(ymin)
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
            mcTotal.Add(obj)
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
    #add legend and LaTeX 
    leg.Draw()
    t1 = rt.TLatex(0.1,0.94, "CMS Preliminary")
    t2 = rt.TLatex(0.55,0.94, "#sqrt{s}=13 TeV, L = "+lumistr)
    t1.SetNDC()
    t2.SetNDC()
    t1.SetTextSize(0.06)
    t2.SetTextSize(0.06)
    t1.Draw()
    t2.Draw()
    #draw data/MC
    c.cd()
    if data and mc: 
        pad2 = rt.TPad("pad2","pad2",0,0.0,1,0.4)
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
