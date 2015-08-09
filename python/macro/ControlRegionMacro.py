import sys
import ROOT as rt

def makeTreeDict(fileDict, treeName, debug=False):
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

def cuts_TTJetsSingleLepton(event, debug=False):
    """Run 2 TTBar single lepton selection"""
    #HLT requirement
    passedTrigger = event.HLTDecision[0] or event.HLTDecision[1] or event.HLTDecision[8] or event.HLTDecision[9]
    if not passedTrigger: return False
    #lepton = tight ele or mu with pt > 30
    if abs(event.lep1Type) != 11 and abs(event.lep1Type) != 13: return False
    if not event.lep1PassTight: return False
    if event.lep1Pt < 30: return False
    #MET and MT cuts
    if event.MET < 30: return False
    if event.lep1MT < 30 or event.lep1MT > 100: return False
    #b-tag requirement
    if event.NBJetsMedium < 1: return False
    #razor baseline cut
    if event.MR < 300 or event.Rsq < 0.15: return False
    #passes selection
    if debug: print("Event passes run 2 TTBarSingleLepton control region selection")
    return True

def fill_RazorBasic(event, hists, weight, debug=False):
    #TODO: check for var in tree
    if "MR" in hists:
        if debug: print("Filling MR histogram")
        hists["MR"].Fill(event.MR, weight)
    if "Rsq" in hists:
        if debug: print("Filling Rsq histogram")
        hists["Rsq"].Fill(event.Rsq, weight)

def standardRun2Weights(event, wHists, debug=False):
    """Apply pileup weights and other known MC correction factors"""
    eventWeight = 1.0
    #TODO: add error checking
    #pileup reweighting
    eventWeight *= wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(event.NPU_0))
    if debug: print("Applying a weight of "+str(eventWeight))
    return eventWeight

def loadRun2WeightHists(debug=False):
    puWeightFileName = "data/Run1PileupWeights.root"
    wHists = {}
    #pileup weight histogram
    if debug: print("Opening pileup weight file "+puWeightFileName)
    puFile = rt.TFile(puWeightFileName)
    wHists["pileup"] = puFile.Get("PUWeight_Run1")
    wHists["pileup"].SetDirectory(0)
    assert wHists["pileup"]
    return wHists

def loopTree(tree, cutf, weightf, weighthists, fillf, hists, maxEvents=-1, debug=False):
    """Loop over a single tree and fill histograms"""
    if debug: print ("Looping tree "+tree.GetName())
    for i,event in enumerate(tree):
        if maxEvents >= 0 and i >= maxEvents: break
        if i % 100000 == 0: print("Processing entry "+str(i))
        if not cutf(event, debug): continue
        w = weightf(event, weighthists, debug)
        fillf(event, hists, w, debug)

def loopTrees(treeDict, cutf, weightf, weighthists, fillf, histDict, maxEvents=-1, debug=False):
    for name in treeDict: 
        loopTree(treeDict[name], cutf, weightf, weighthists, fillf, histDict[name], maxEvents, debug)

def makeStack(hists, ordering, title="Stack"):
    stack = rt.THStack("thstack"+title.replace(" ",""), title)
    for name in ordering: 
        stack.Add(hists[name])
        #TODO: error checking
    return stack

def makeLegend(hists, names, ordering, x1=0.6, y1=0.6, x2=0.9, y2=0.9):
    leg = rt.TLegend(x1, y1, x2, y2)
    for name in ordering: 
        leg.AddEntry(hists[name], names[name])
    return leg

def basicplot(mc=0, data=0, fit=0, leg=0, xtitle="", ytitle="Number of events", ymin=0.1, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", ratiomin=0.5, ratiomax=1.5, saveroot=False, savepdf=False, savepng=True):
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
        mc.GetYaxis().SetTitleOffset(0.45)
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

def plotmacro_test(histDict):
    ordering = ["SingleTop"]

    names = {}
    MRHists = {}
    for name in histDict: names[name] = name
    for name in histDict: MRHists[name] = histDict[name]["MR"]

    stack = makeStack(MRHists, ordering, "MR")
    legend = makeLegend(MRHists, names, reversed(ordering))
    basicplot(mc=stack, leg=legend, xtitle="MR", printstr="MR_TTJetsSingleLepton")

def sethistcolor(hist, name):
    if name == "SingleTop": hist.SetFillColor(rt.kOrange)
    else: hist.SetFillColor(kBlack)

if __name__ == "__main__":
    rt.gROOT.SetBatch()
    if len(sys.argv) > 1: debug = bool(sys.argv[1])
    else: debug = False
    maxEvents = 10000

    #setup
    inputs = {"SingleTop":"SingleTop_1pb_weighted.root"}
    #load files
    files = {}
    for name in inputs: 
        files[name] = rt.TFile(inputs[name])
        assert files[name]
    #load trees
    trees = makeTreeDict(files, "ControlSampleEvent", debug)
    #load weight histograms
    weighthists = loadRun2WeightHists(debug)

    #define histograms to fill
    hists = {"SingleTop": {}}
    for name in hists:
        hists[name]["MR"] = rt.TH1F("MR"+name, "M_{R} (GeV)", 20, 300, 4000)
        for var in hists[name]: 
            hists[name][var].Sumw2()
            sethistcolor(hists[name][var], name)
    
    loopTrees(trees, cuts_TTJetsSingleLepton, standardRun2Weights, weighthists, fill_RazorBasic, hists, maxEvents=maxEvents, debug=debug)
    plotmacro_test(hists)

#TODO: standardize naming convention
#TODO: add debugging info to each function
#TODO: split into multiple files
#TODO: add other TTJets samples
#TODO: expand to other CRs
#TODO: add fit
#TODO: comment everything
