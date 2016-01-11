import os
import ROOT as rt
import copy
from array import *

#local imports
from PlotFit import setFFColors
from RunCombine import exec_me
from razorWeights import applyMTUncertainty1D, applyMTUncertainty2D

def makeLegend(hists, titles, ordering, x1=0.6, y1=0.6, x2=0.9, y2=0.9):
    """Takes a dict of histograms, a dict of histogram titles, and an ordered list of names, and returns a legend with the histograms in the desired order"""
    leg = rt.TLegend(x1, y1, x2, y2)
    rt.SetOwnership(leg, False)
    for name in ordering: 
        leg.AddEntry(hists[name], titles[name])
    return leg

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

def getLinesForUnrolled(hist):        
    # the gray lines
    lines = []
    nBinsX = hist.GetNbinsX()
    nBinsY = hist.GetNbinsY()
    for i in range(1,nBinsX):
        lineX = i*nBinsY
        lines.append(rt.TLine(lineX, 0, lineX, hist.GetYaxis().GetXmax()))
        lines[-1].SetLineStyle(2)
        lines[-1].SetLineColor(rt.kBlack)
    return lines

def setHistColor(hist, name):
    colors = {"WJets":rt.kRed+1, "WJetsInv":rt.kRed+1, "DYJets":rt.kBlue+1, "DYJetsInv":rt.kBlue+1, "TTJets":rt.kGreen+2, "ZInv":rt.kCyan+1, "QCD":rt.kMagenta, "SingleTop":rt.kOrange-3, "VV":rt.kViolet+3, "TTV":rt.kGreen-7, "DYJetsLow":rt.kBlue+1, "PhotonInv":rt.kOrange , "Other":rt.kAzure+4}
    """Sets histogram color"""
    if name in colors: hist.SetFillColor(colors[name])
    else: print "Warning in macro.py: histogram fill color not set for",name


def plot_several(c, hists=0, leg=0, colors=[], xtitle="", ytitle="Events", ymin=None, ymax=None, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", commentstr="", ratiomin=0.5, ratiomax=1.5, saveroot=True, savepdf=True, savepng=True, printdir='.'):
    """Draw several histograms as colored lines"""
    #setup
    c.Clear()
    c.cd()
    if len(hists) < 1: 
        print "Error in plot_several: no histograms provided!"
        return
    elif len(hists) > 1: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.4, 1, 1)
    else: pad1 = rt.TPad(printstr+"pad1", printstr+"pad1", 0, 0.1, 1, 1)
    pad1.SetBottomMargin(0)
    pad1.SetLogx(logx)
    pad1.SetLogy(logy)
    pad1.Draw()
    pad1.cd()

    #draw 
    for i in range(len(hists)-1,-1,-1): #draw in reverse order
        hist = hists[i]
        if i < len(colors):
            hist.SetMarkerColor(colors[i])
        hist.SetLineWidth(2)
        hist.SetTitle("")
        hist.GetYaxis().SetTitle(ytitle)
        hist.GetYaxis().SetLabelSize(0.03)
        hist.SetStats(0)
        if logy: hist.GetXaxis().SetMoreLogLabels()
        if len(hists) > 1:
            hist.GetYaxis().SetTitleOffset(0.50)
        else:
            hist.GetXaxis().SetTitle(xtitle)
            hist.GetYaxis().SetTitleOffset(0.45)
        hist.GetYaxis().SetTitleSize(0.05)
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(0.5)
        if ymin is not None: hist.SetMinimum(ymin)
        if ymax is not None: hist.SetMaximum(ymax)
        if i == len(hists)-1:
            hist.Draw('p')
        else:
            hist.Draw("psame")
    pad1.Modified()
    rt.gPad.Update()

    #add legend and LaTeX 
    leg.Draw()
    t1 = rt.TLatex(0.1,0.94, "CMS preliminary")
    t2 = rt.TLatex(0.55,0.94, ((lumistr != "")*((lumistr)+' ('))+'13 TeV'+((lumistr != '')*(')')))
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

    #ratio histograms
    c.cd()
    pad2 = rt.TPad(printstr+"pad2",printstr+"pad2",0,0.0,1,0.4)
    pad2.SetTopMargin(0)
    pad2.SetTopMargin(0.008)
    pad2.SetBottomMargin(0.25)
    pad2.SetGridy()
    pad2.SetLogx(logx)
    lowerHists = []
    for i in range(len(hists)-1,-1,-1): #draw in reverse order
        if i == 0: continue
        hist = hists[i]
        lowerHists.append(make1DRatioHistogram(hist, hists[0], xtitle, ratiomin, ratiomax, logx))
        lowerHists[-1].SetMarkerSize(0.5)
        lowerHists[-1].SetMarkerStyle(20)
        lowerHists[-1].GetYaxis().SetTitle("Ratio")

        if i == len(hists)-1:
            pad2.Draw()
            pad2.cd()
            lowerHists[-1].Draw('p') 
        else:
            lowerHists[-1].Draw('psame')
        pad2.Modified()
        rt.gPad.Update()

    #save
    if savepng: c.Print(printdir+'/'+printstr+".png")
    if savepdf: c.Print(printdir+'/'+printstr+".pdf")
    if saveroot: c.Print(printdir+'/'+printstr+".root")

def plot_basic(c, mc=0, data=0, fit=0, leg=0, xtitle="", ytitle="Events", ymin=None, ymax=None, printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", commentstr="", ratiomin=0.5, ratiomax=1.5, pad2Opt="Ratio", fitColor=rt.kRed, mcErrColor=rt.kBlack, customPad2Hist=None, saveroot=True, savepdf=True, savepng=True, printdir='.', grayLines=None):
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
        mc.SetTitle("")
        #make total MC histogram
        histList = mc.GetHists()
        mcTotal = histList.First().Clone()
        mcTotal.SetTitle("")
        mcTotal.SetStats(0)
        mcTotal.Reset()
        for h in histList:
            mcTotal.Add(h)
        mcTotal.SetFillColor(mcErrColor)
        mcTotal.GetYaxis().SetTitle(ytitle)
        mcTotal.GetYaxis().SetTitleSize(0.06)
        mcTotal.GetYaxis().SetLabelSize(0.06)
        if logy: mcTotal.GetXaxis().SetMoreLogLabels()
        if not data: mcTotal.GetXaxis().SetTitle(xtitle)
        mcTotal.GetYaxis().SetTitleOffset(0.6)
        mcTotal.GetYaxis().SetTitleSize(0.05)
        if ymin is not None: mcTotal.SetMinimum(ymin)
        if ymax is not None: mcTotal.SetMaximum(ymax)
        if data and data.GetMaximum() > mcTotal.GetMaximum() and ymax is None: 
            mcTotal.SetMaximum(data.GetMaximum())
        mc.Draw("hist")
        mc.GetXaxis().SetTitle(xtitle)
        mc.GetYaxis().SetTitle(ytitle)
        mc.GetYaxis().SetTitleOffset(0.60)
        mc.GetYaxis().SetTitleSize(0.06)
        mc.GetYaxis().SetLabelSize(0.06)
        mcTotal.SetFillStyle(3001)
        mcTotal.Draw("e2same")
    #draw fit
    if fit:
        fit.SetStats(0)
        fit.SetMarkerStyle(21)
        fit.SetLineWidth(2)
        if mc:
            fit.SetMarkerSize(1.2)
            fit.SetLineColor(fitColor+2)
        else:
            fit.SetMarkerSize(0)
            fit.SetLineColor(fitColor)
        fit.SetMarkerColor(fitColor)
        fit.SetFillStyle(3001)
        fit.SetFillColor(rt.kAzure+6)
        fit.GetYaxis().SetTitle(ytitle)
        fit.SetTitle("")
        if ymin is not None: fit.SetMinimum(ymin)
        if ymax is not None: fit.SetMaximum(ymax)
        if mc:
            fit.Draw("pesame")
        else:
            fit.Draw("e2same")
            fitCopy = fit.Clone()
            fitCopy.SetFillStyle(0)        
            fitCopy.SetLineWidth(2)
            fitCopy.Draw("histsame")
    #draw data
    if data:
        data.SetStats(0)
        data.SetMarkerStyle(20)
        data.SetMarkerSize(1)
        data.SetLineWidth(2)
        data.SetLineColor(rt.kBlack)
        data.GetYaxis().SetTitle(ytitle)
        data.SetTitle("")
        if ymin is not None: data.SetMinimum(ymin)
        if ymax is not None: data.SetMaximum(ymax)
        data.SetBinErrorOption(rt.TH1.kPoisson)
        data.Draw("pe1same")
    pad1.Modified()
    rt.gPad.Update()

    #gray lines
    if grayLines is not None:
        for line in grayLines: 
            line.Draw()

    #add legend and LaTeX 
    leg.Draw()
    t1 = rt.TLatex(0.1,0.94, "CMS preliminary")
    t2 = rt.TLatex(0.55,0.94, ((lumistr != "")*(lumistr)+' (')+'13 TeV'+((lumistr != '')*(')')))
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
    lowerPadHist2 = None

    #set up custom histogram if provided
    if customPad2Hist is not None:
        customPad2Hist.SetTitle("")
        customPad2Hist.GetXaxis().SetTitle(xtitle)
        customPad2Hist.SetMinimum(ratiomin)
        customPad2Hist.SetMaximum(ratiomax)
        customPad2Hist.SetStats(0)
        customPad2Hist.GetXaxis().SetLabelSize(0.06)
        customPad2Hist.GetYaxis().SetLabelSize(0.06)
        customPad2Hist.GetYaxis().SetTitleOffset(0.35)
        customPad2Hist.GetXaxis().SetTitleOffset(1.50)
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
            #draw data/MC with poisson errors from data
            lowerPadHist = make1DRatioHistogram(data, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreDenominatorErrs=True)
            lowerPadHist.GetYaxis().SetTitle("Data / MC")
            #draw MC uncertainties on the ratio
            #get plot style from MC histogram
            dataWithMCStyle = mcTotal.Clone()
            for bx in range(1, dataWithMCStyle.GetNbinsX()+1):
                dataWithMCStyle.SetBinContent(bx, data.GetBinContent(bx))
                dataWithMCStyle.SetBinError(bx, data.GetBinError(bx))
            lowerPadHist2 = make1DRatioHistogram(dataWithMCStyle, mcTotal, xtitle, ratiomin, ratiomax, logx, ignoreNumeratorErrs=True)
            lowerPadHist2.GetYaxis().SetTitle(lowerPadHist.GetYaxis().GetTitle())
            for bx in range(1, lowerPadHist2.GetNbinsX()+1):
                lowerPadHist2.SetBinContent(bx,1)
            lowerPadHist2.SetFillStyle(3001)
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
            lowerPadHist.GetYaxis().SetTitle("Nsigma (Stat+Sys)")
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
    if lowerPadHist2 is not None:
        pad2.Draw()
        pad2.cd()
        lowerPadHist2.Draw("e2same")
        pad2.Modified()
        rt.gPad.Update()
    if lowerPadHist is not None:
        pad2.Draw()
        pad2.cd()
        lowerPadHist.Draw("pesame")
        pad2.Modified()
        rt.gPad.Update()

    #save
    if savepng: c.Print(printdir+'/'+printstr+".png")
    if savepdf: c.Print(printdir+'/'+printstr+".pdf")
    if saveroot: c.Print(printdir+'/'+printstr+".root")


def draw2DHist(c, hist, xtitle="", ytitle="", ztitle="", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, drawErrs=False, palette=53, grayGraphs=None, saveroot=True, savepdf=True, savepng=True, numDigits=1, printdir='.'):
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
    t1 = rt.TLatex(0.1,0.94, "CMS preliminary")
    t2 = rt.TLatex(0.55,0.94, ((lumistr != "")*(lumistr)+' (')+'13 TeV'+((lumistr != '')*(')')))
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

def make1DRatioHistogram(num, denom, xtitle="", ratiomin=0.25, ratiomax=2.0, logx=False, forPad2=True, ignoreNumeratorErrs=False, ignoreDenominatorErrs=False):
    ratio = num.Clone()
    denomClone = denom.Clone()
    if ignoreNumeratorErrs:
        for bx in range(1, ratio.GetNbinsX()+1):
            ratio.SetBinError(bx, 0.0)
    if ignoreDenominatorErrs:
        for bx in range(1, denomClone.GetNbinsX()+1):
            denomClone.SetBinError(bx, 0.0)
    ratio.Divide(denomClone)
    ratio.SetTitle("")
    ratio.GetXaxis().SetTitle(xtitle)
    ratio.SetMinimum(ratiomin)
    ratio.SetMaximum(ratiomax)
    ratio.SetStats(0)
    if forPad2:
        ratio.GetXaxis().SetLabelSize(0.1)
        ratio.GetYaxis().SetLabelSize(0.08)
        ratio.GetYaxis().SetTitleOffset(0.35)
        ratio.GetXaxis().SetTitleOffset(1.50)
        ratio.GetYaxis().SetTitleSize(0.08)
        ratio.GetXaxis().SetTitleSize(0.08)
        ratio.SetTitle("")
        #if logx: ratio.GetXaxis().SetMoreLogLabels()
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
        ret.GetXaxis().SetTitleOffset(1.50)
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
        ret.GetXaxis().SetTitleOffset(1.50)
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

def plot_basic_2D(c, mc=0, data=0, fit=0, xtitle="", ytitle="", ztitle="Events", zmin=None, zmax=None, printstr="hist", logx=True, logy=True, logz=True, lumistr="", commentstr="", dotext=True, saveroot=True, savepdf=True, savepng=True, nsigmaFitData=None, nsigmaFitMC=None, mcDict=None, mcSamples=None, printdir="."):
    """Plotting macro for data, MC, and/or fit yields.  Creates french flag plots comparing data/MC/fit if able.
    mc, data, fit: 2d histograms for MC prediction, data, and fit prediction (all are optional)
    lumistr: tex-formatted string indicating the integrated luminosity
    commentstr: additional string that will be displayed at the top of the plot 
    nsigmaFitData, nsigmaFitMC (optional): externally provided nsigma histograms to use instead of (yield-prediction)/sigma 
    mcDict (optional): dictionary of MC histograms, for making stacked unrolled plots (provide a list of sample names, mcSamples, to enforce an ordering on the MC histograms in the stack)
    """
    #make a gray square for each -999 bin
    grayGraphs = [makeGrayGraphs(hist) for hist in [mc,fit,data]]
    #unroll 2D hists to plot in 1D
    unrolled = unroll2DHistograms([mc, data, fit])
    #if individual MC histograms are available, unroll all of them
    if mcDict is not None:
        if mcSamples is not None:
            mcUnrolledList = unroll2DHistograms([mcDict[s] for s in mcSamples])
            mcUnrolledDict = {mcSamples[n]:mcUnrolledList[n] for n in range(len(mcSamples))}
        else:
            mcSamples = [s for s in mcDict]
            mcUnrolledList = unroll2DHistograms([mcDict[s] for s in mcSamples])
            mcUnrolledDict = {mcSamples[n]:mcUnrolledList[n] for n in range(len(mcSamples))}
        for s in mcSamples: setHistColor(mcUnrolledDict[s], s)
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
        if mcDict is not None: #draw full stack of MC on unrolled plot
            mcStack = makeStack(mcUnrolledDict, mcSamples, "MC")
        else: #draw sum of MC only
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
            #(first remove any blinded bins)
            legDataMC = rt.TLegend(0.7, 0.7, 0.9, 0.9)
            rt.SetOwnership(legDataMC, False)
            if mcDict is not None:
                blindMCUnrolledDict = {}
                #blind and create stack
                for s in mcSamples: 
                    blindMCUnrolledDict[s] = mcUnrolledDict[s].Clone()
                    for bx in range(1,blindMCUnrolledDict[s].GetNbinsX()+1):
                        if unrolled[1].GetBinContent(bx) < 0:
                            blindMCUnrolledDict[s].SetBinContent(bx, -999)
                blindStack = makeStack(blindMCUnrolledDict, mcSamples, "MC")
                #add each sample to legend
                for n in range(len(mcSamples)-1,-1,-1):
                    s = mcSamples[n]
                    legDataMC.AddEntry(blindMCUnrolledDict[s], s)
            else:
                blindMC = unrolled[0].Clone("blindMC")
                for bx in range(1,blindMC.GetNbinsX()+1):
                    if unrolled[1].GetBinContent(bx) < 0:
                        blindMC.SetBinContent(bx,-999)
                blindStack = makeStack({"MC":blindMC}, ["MC"], "MC")
                legDataMC.AddEntry(blindMC, "MC Prediction")
            legDataMC.AddEntry(unrolled[1], "Data")
            plot_basic(c, blindStack, unrolled[1], None, legDataMC, xtitle="Bin", ymin=0.1, printstr=printstr+"UnrolledDataMC", lumistr=lumistr, mcErrColor=rt.kBlack, commentstr=commentstr, ratiomin=0.5,ratiomax=1.5, pad2Opt="ratio", saveroot=True, printdir=printdir)
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
            nsigmaUnrolledFitMC = None
            if nsigmaFitMC is not None:
                nsigmaUnrolledFitMC = unroll2DHistograms([nsigmaFitMC])[0]
                for bx in range(1, nsigmaUnrolledFitMC.GetNbinsX()+1):
                    nsigmaUnrolledFitMC.SetBinError(bx,0.0)

            legMCFit = rt.TLegend(0.7, 0.7, 0.9, 0.9)
            rt.SetOwnership(legMCFit, False)
            legMCFit.AddEntry(unrolled[0], "MC Prediction")
            legMCFit.AddEntry(unrolled[2], "Fit Prediction")
            plot_basic(c, mcStack, None, unrolled[2], legMCFit, xtitle="Bin", ymin=0.1, printstr=printstr+"UnrolledMCFit", lumistr=lumistr, commentstr=commentstr, ratiomin=-5., ratiomax=5.0, pad2Opt="ff", saveroot=True, mcErrColor=rt.kBlack, customPad2Hist=nsigmaUnrolledFitMC, printdir=printdir)
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

            #unroll and compare
            grayLines = getLinesForUnrolled(fit)

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

            plot_basic(c, None, unrolled[1], blindFit, legDataFit, xtitle="Bin", ymin=0.1, printstr=printstr+"UnrolledDataFit", lumistr=lumistr, commentstr=commentstr, ratiomin=-5., ratiomax=5.0, pad2Opt="ff", fitColor=rt.kBlue, saveroot=True, customPad2Hist=nsigmaUnrolled, printdir=printdir, grayLines=grayLines)
    if fit:
        draw2DHist(c, fit, xtitle, ytitle, ztitle, zmin, zmax, printstr+'Fit', lumistr=lumistr, commentstr=commentstr+", Fit prediction", grayGraphs=grayGraphs[1], dotext=dotext, drawErrs=True, saveroot=saveroot, savepdf=savepdf, savepng=savepng, printdir=printdir)
        #make relative uncertainty histogram
        relUncFitHist = make2DRelativeUncertaintyHistogram(fit)
        draw2DHist(c, relUncFitHist, xtitle, ytitle, ztitle, 0.0, 2.0, printstr+'FitRelUnc', lumistr=lumistr, commentstr='#splitline{'+commentstr+"}{#it{fit relative uncertainty}}", grayGraphs=grayGraphs[1], dotext=dotext, drawErrs=False, saveroot=saveroot, savepdf=savepdf, savepng=savepng, logz=False, printdir=printdir)

def makeStackAndPlot(canvas, mcHists={}, dataHist=None, dataName="Data", mcOrdering=[], titles=[], mcTitle="Stack", xtitle="", ytitle="Events", printstr="hist", logx=False, logy=True, lumistr="40 pb^{-1}", saveroot=True, savepdf=True, savepng=True, ymin=None, ymax=None):
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

def table_basic(headers=[], cols=[], caption="", printstr='table', landscape=False, printdir='.'):
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

    with open(printdir+'/'+printstr+'.tex', 'w') as f:
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

