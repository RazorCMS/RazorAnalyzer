from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
from WriteDataCard import *
import os
import random
import sys
import math


def setStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetPaintTextFormat("1.2g")

def findLastBin(h):
    for i in range(1,h.GetXaxis().GetNbins()):
        thisbin = h.GetXaxis().GetNbins()-i
        if h.GetBinContent(thisbin)>=0.1: return thisbin+1
    return h.GetXaxis().GetNbins()

def getPads(c):    
    pad1 = rt.TPad(c.GetName()+"_pad1","pad1",0,0.25,1,1)
    pad2 = rt.TPad(c.GetName()+"_pad2","pad2",0,0,1,0.25)
    pad1.Range(-213.4588,-0.3237935,4222.803,5.412602);
    pad2.Range(-213.4588,-2.206896,4222.803,3.241379);
    pad1.SetLeftMargin(0.15)
    pad2.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad2.SetRightMargin(0.05)
    pad1.SetTopMargin(0.12)
    pad2.SetTopMargin(0.)
    pad1.SetBottomMargin(0.)
    pad2.SetBottomMargin(0.47)
    pad1.Draw()
    pad1.cd()
    rt.gPad.SetLogy(1)
    return pad1, pad2


def setDataHist(h_data,xTitle,yTitle,color=rt.kBlack):        
    h_data.SetMarkerColor(color)
    h_data.SetMarkerStyle(20)
    h_data.SetLineColor(color)
    h_data.GetXaxis().SetTitle(xTitle)
    h_data.GetYaxis().SetTitle(yTitle)
    h_data.GetXaxis().SetLabelOffset(0.16)
    h_data.GetXaxis().SetLabelSize(0.06)
    h_data.GetYaxis().SetLabelSize(0.06)
    h_data.GetXaxis().SetTitleSize(0.06)
    h_data.GetYaxis().SetTitleSize(0.08)
    h_data.GetXaxis().SetTitleOffset(0.8)
    h_data.GetYaxis().SetTitleOffset(0.7)
    h_data.GetXaxis().SetTicks("+-")
    h_data.SetMaximum(math.pow(h_data.GetBinContent(h_data.GetMaximumBin()),1.25))
    h_data.SetMinimum(max(1e-1,1e-1*h_data.GetBinContent(h_data.GetMinimumBin())))
    return h_data

def getDivideHistos(h,hClone,h_data,xTitle,divTitle):
    h.Sumw2()
    hClone.Sumw2()
    h_data.Sumw2()

    hDivide = h.Clone(h.GetName()+"Divide") 
    hCloneDivide = hClone.Clone(hClone.GetName()+"Divide") 
    hDataDivide = h_data.Clone(h_data.GetName()+"Divide")
    hDivide.Sumw2()
    hCloneDivide.Sumw2()
    hDataDivide.Sumw2()
    hCloneDivide.GetYaxis().SetLabelSize(0.18)
    hCloneDivide.SetTitle("")
    hCloneDivide.SetMaximum(3.5)
    hCloneDivide.SetMinimum(0.)
    hCloneDivide.GetXaxis().SetLabelSize(0.22)
    hCloneDivide.GetXaxis().SetTitleSize(0.22)

    
    for i in range(1, h_data.GetNbinsX()+1):
        tmpVal = hCloneDivide.GetBinContent(i)
        if tmpVal != -0.:
            hDataDivide.SetBinContent(i, hDataDivide.GetBinContent(i)/tmpVal)
            hDataDivide.SetBinError(i, hDataDivide.GetBinError(i)/tmpVal)
            hCloneDivide.SetBinContent(i, hCloneDivide.GetBinContent(i)/tmpVal)
            hCloneDivide.SetBinError(i, hCloneDivide.GetBinError(i)/tmpVal)
            hDivide.SetBinContent(i, hDivide.GetBinContent(i)/tmpVal)
            hDivide.SetBinError(i, hDivide.GetBinError(i)/tmpVal)

            
    hCloneDivide.GetXaxis().SetTitleOffset(0.97)
    hCloneDivide.GetXaxis().SetLabelOffset(0.02)
    hCloneDivide.GetXaxis().SetTitle(xTitle)

    hCloneDivide.GetYaxis().SetNdivisions(504,rt.kTRUE)
    hCloneDivide.GetYaxis().SetTitleOffset(0.2)
    hCloneDivide.GetYaxis().SetTitleSize(0.22)
    hCloneDivide.GetYaxis().SetTitle(divTitle)
    hCloneDivide.GetXaxis().SetTicks("+")
    hCloneDivide.GetXaxis().SetTickLength(0.07)
    hCloneDivide.SetMarkerColor(hCloneDivide.GetFillColor())
    
    return hDivide, hCloneDivide, hDataDivide
    
def print1DProj(c,h,h_data,printName,xTitle,yTitle,lumiLabel="",boxLabel="",tLeg=None,h_components=[],h_colors=[],h_labels=[]):
    
    pad1, pad2 = getPads(c)

    h.SetLineWidth(2)
    h.SetLineColor(rt.kBlue)
    hClone = h.Clone(h.GetName()+"Clone")
    hClone.SetLineColor(rt.kBlue)
    hClone.SetFillColor(rt.kBlue-10)
    
    h_data = setDataHist(h_data,xTitle,yTitle)
    
    h_data.Draw("pe")
    hClone.Draw("e2same")
    h.SetFillStyle(0)
    for h_comp, color, label in zip(h_components, h_colors, h_labels):
        h_comp.SetLineColor(color)
        h_comp.SetLineWidth(2)            
        h_comp.Draw("histsame")
    h.DrawCopy("histsame")
    h_data.Draw("pesame")
    pad1.Draw()
    c.Update()
    c.cd()
    pad2.Draw()
    pad2.cd()
    rt.gPad.SetLogy(0)

    hDivide, hCloneDivide, hDataDivide  = getDivideHistos(h, hClone, h_data, xTitle, "Data/Fit")
    
    hCloneDivide.Draw("e2")
    #hDivide.Draw("histsame")
    hDataDivide.Draw('pesame')
    hCloneDivide.Draw("axissame")


    pad2.Update()
    pad1.cd()
    pad1.Update()
    pad1.Draw()

    if tLeg==None:
        if len(h_components)>=7:
            tLeg = rt.TLegend(0.7,0.3,0.9,0.8)
        elif len(h_components)==6:
            tLeg = rt.TLegend(0.7,0.35,0.9,0.8)
        elif len(h_components)==5:
            tLeg = rt.TLegend(0.7,0.4,0.9,0.8)
        elif len(h_components)==4:
            tLeg = rt.TLegend(0.7,0.45,0.9,0.8)
        elif len(h_components)==3:
            tLeg = rt.TLegend(0.7,0.5,0.9,0.8)
        elif len(h_components)==2:
            tLeg = rt.TLegend(0.7,0.55,0.9,0.8)            
        else:
            tLeg = rt.TLegend(0.7,0.6,0.9,0.8)
        tLeg.SetFillColor(0)
        tLeg.SetTextFont(42)
        tLeg.SetLineColor(0)
        tLeg.AddEntry(h_data,"Sim. Data","lep")
        tLeg.AddEntry(hClone,"Fit Total","lf")
        for h_comp, color, label in zip(h_components, h_colors, h_labels):                
            tLeg.AddEntry(h_comp,label,"l")
            
    tLeg.Draw("same")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.9,"CMS simulation")
    l.DrawLatex(0.78,0.9,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.DrawLatex(0.2,0.8,boxLabel)

    c.cd()
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')


def print1DSlice(c,h_slices,h_data_slices,printName,xTitle,yTitle,lumiLabel="",boxLabel="",tLeg=None,h_colors=[],h_labels=[]):

    pad1, pad2 = getPads(c)

    for h,color in zip(h_slices,h_colors):
        h.SetLineWidth(2)
        h.SetLineColor(color)

    for h_data,color in zip(h_data_slices,h_colors):
        h_data = setDataHist(h_data,xTitle,yTitle,color)
    #h_data_slices[0].SetMaximum(1e2*h_data_slices[0].GetMaximum())
            
    first = True
    for h_data,color in zip(h_data_slices,h_colors):
        if first:            
            h_data.Draw("pe")
            first = False
        else:
            h_data.Draw("pesame")
    
    for h,color in zip(h_slices,h_colors): h.Draw("histsame")
        
    for h_data, color in zip(h_data_slices,h_colors): h_data.Draw("pesame")
        
    pad1.Draw()
    c.Update()
    c.cd()
    pad2.Draw()
    pad2.cd()
    rt.gPad.SetLogy(0)

    hDataDivides=[]
    hDivides = []
    hCloneDivides = []
    for h, h_data, color in zip(h_slices, h_data_slices, h_colors):
        hClone = h.Clone(h.GetName()+"Clone")
        hDivide, hCloneDivide, hDataDivide  = getDivideHistos(h, hClone, h_data, xTitle, "Data/Fit")
        hCloneDivide.SetLineColor(rt.kWhite)
        
        hDataDivides.append(hDataDivide)
        hCloneDivides.append(hCloneDivide)
        hDivides.append(hDivide)

    hCloneDivides[0].Draw("axis")
    hDataDivides[0].Draw("pesame")
    for i, hDataDivide in enumerate(hDataDivides):
        if i>0: hDataDivide.Draw("pesame")
    #hCloneDivides[0].Draw("axissame")
    
    pad2.Update()
    pad1.cd()
    pad1.Update()
    pad1.Draw()

    if tLeg==None:
        if len(h_slices)>=7:
            tLeg = rt.TLegend(0.7,0.3,0.9,0.8)
        elif len(h_slices)==6:
            tLeg = rt.TLegend(0.7,0.35,0.9,0.8)
        elif len(h_slices)==5:
            tLeg = rt.TLegend(0.7,0.4,0.9,0.8)
        elif len(h_slices)==4:
            tLeg = rt.TLegend(0.7,0.45,0.9,0.8)
        elif len(h_slices)==3:
            tLeg = rt.TLegend(0.7,0.5,0.9,0.8)
        elif len(h_slices)==2:
            tLeg = rt.TLegend(0.7,0.55,0.9,0.8)            
        else:
            tLeg = rt.TLegend(0.7,0.6,0.9,0.8)
        tLeg.SetFillColor(0)
        tLeg.SetTextFont(42)
        tLeg.SetLineColor(0)
        for h, color, label in zip(h_slices, h_colors, h_labels):                
            tLeg.AddEntry(h,label,"l")
            
    tLeg.Draw("same")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.9,"CMS simulation")
    l.DrawLatex(0.78,0.9,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.DrawLatex(0.2,0.8,boxLabel)

    c.cd()
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    
def setFFColors(hNS, minZ=-5.1, maxZ=5.1):
    Red = array('d',  [0.00, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00])
    Green = array('d',[0.00, 0.70, 0.90, 1.00, 0.90, 0.70, 0.00])
    Blue = array('d', [1.00, 1.00, 1.00, 1.00, 0.90, 0.70, 0.00])
    Length =array('d',[0.00, 0.20, 0.35, 0.50, 0.65, 0.8, 1.00]) # colors get darker faster at 4sigma
    rt.TColor.CreateGradientColorTable(7,Length,Red,Green,Blue,999)
    hNS.SetMaximum(maxZ)
    hNS.SetMinimum(minZ) # so the binning is 0 2 4
    hNS.SetContour(999)

    
def setRainbowColors(hNS, minZ=0, maxZ=100):
    rt.gStyle.SetPalette(1)
    Red = array('d',  [0.00, 0.00, 0.87, 1.00, 0.51])
    Green = array('d',[0.00, 0.81, 1.00, 0.20, 0.00])
    Blue = array('d', [0.51, 1.00, 0.12, 0.00, 0.00])
    Length =array('d',[0.00, 0.34, 0.61, 0.84, 1.00]) 
    rt.TColor.CreateGradientColorTable(5,Length,Red,Green,Blue,999)
    hNS.SetMaximum(maxZ)
    hNS.SetMinimum(minZ) # so the binning is 0 2 4
    hNS.SetContour(999)

    
def print2DCanvas(c,h,printName):
    c.SetLogx(1)
    c.SetLogy(1)
    h.Draw("textcolz")
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')

def set2DHisto(h2D,xTitle,yTitle,zTitle):    
    h2D.GetXaxis().SetMoreLogLabels()
    h2D.GetXaxis().SetNoExponent()
    h2D.GetYaxis().SetMoreLogLabels()
    h2D.GetYaxis().SetNoExponent()
    h2D.SetMarkerSize(1.5)
    h2D.GetXaxis().SetTitle(xTitle)
    h2D.GetYaxis().SetTitle(yTitle)
    h2D.GetZaxis().SetTitle(zTitle)    
    h2D.GetXaxis().SetTitleSize(0.056)
    h2D.GetXaxis().SetLabelSize(0.056)
    h2D.GetYaxis().SetTitleSize(0.056)
    h2D.GetYaxis().SetLabelSize(0.056)
    h2D.GetZaxis().SetLabelSize(0.056)
    h2D.GetZaxis().SetTitleSize(0.056)
    h2D.GetZaxis().SetTitleOffset(1)
    #h2D.GetXaxis().SetTitleOffset(0.8)
    #h2D.GetYaxis().SetTitleOffset(0.7)
    #h2D.GetXaxis().SetTicks("+-")
    return h2D

def getGrayLines(x,y):        
    # the gray lines
    xLines = []
    yLines = []

    lastX = len(x)-1
    lastY = len(y)-1

    for i in range(1,lastY):
        xLines.append(rt.TLine(x[0], y[i], x[lastX], y[i]))
        xLines[i-1].SetLineStyle(2);
        xLines[i-1].SetLineColor(rt.kGray);
        
    for i in range(1,lastX):
        yLines.append(rt.TLine(x[i], y[0], x[i], y[lastY]))
        yLines[i-1].SetLineStyle(2)
        yLines[i-1].SetLineColor(rt.kGray)
        
    return xLines,yLines

def set2DCanvas(c):
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    c.SetRightMargin(0.17)
    c.SetLogx(1)
    c.SetLogy(1)
    c.SetLogz(0)
    return c

def print2DScatter(c,h,printName,xTitle,yTitle,zTitle,lumiLabel,boxLabel,x,y,zMin,zMax):

    c = set2DCanvas(c)
    c.SetLogz(1)
    
    h = set2DHisto(h,xTitle,yTitle,zTitle)
    setRainbowColors(h,zMin,zMax)
    
    h.Draw("colz")
    xLines, yLines = getGrayLines(x,y)
    
    [xLine.Draw("l") for xLine in xLines]
    [yLine.Draw("l") for yLine in yLines]

    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.91,"CMS simulation")
    l.DrawLatex(0.65,0.91,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.DrawLatex(0.2,0.85,boxLabel)
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')

def print2DResiduals(c,h,h_data,printName,xTitle,yTitle,zTitle,lumiLabel,boxLabel,x,y):

    c = set2DCanvas(c)
    h_resi = h_data.Clone(h_data.GetName()+"residuals")
    h_resi.Add(h,-1.)
    h_resi = set2DHisto(h_resi,xTitle,yTitle,zTitle)    
    absMax = max(abs(h_resi.GetMinimum()),abs(h_resi.GetMaximum()))
    setFFColors(h_resi,-1.5*absMax,1.5*absMax)
    h_resi.Draw("colztext")
    xLines, yLines = getGrayLines(x,y)
    
    [xLine.Draw("l") for xLine in xLines]
    [yLine.Draw("l") for yLine in yLines]
        
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.91,"CMS simulation")
    l.DrawLatex(0.65,0.91,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.DrawLatex(0.2,0.85,boxLabel)
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    
def print2DPercentDiff(c,h,h_data,printName,xTitle,yTitle,zTitle,lumiLabel,boxLabel,x,y):
    
    c = set2DCanvas(c)
    h_resi = h_data.Clone(h_data.GetName()+"percentdiff")
    h_resi.Add(h,-1.)
    h_resi.Divide(h)
    h_resi = set2DHisto(h_resi,xTitle,yTitle,zTitle)
    absMax = max(abs(h_resi.GetMinimum()),abs(h_resi.GetMaximum()))
    setFFColors(h_resi,-1.5*absMax,1.5*absMax)
    h_resi.Draw("colztext")
    xLines, yLines = getGrayLines(x,y)
    
    [xLine.Draw("l") for xLine in xLines]
    [yLine.Draw("l") for yLine in yLines]

    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.91,"CMS simulation")
    l.DrawLatex(0.65,0.91,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.DrawLatex(0.2,0.85,boxLabel)
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    
def print2DNSigma(c,h,h_data,printName,xTitle,yTitle,zTitle,lumiLabel,boxLabel,x,y):
    
    c = set2DCanvas(c)
    h_resi = h_data.Clone(h_data.GetName()+"nsigma")
    h_resi.Add(h,-1.)
    for i in range(1,h_resi.GetNbinsX()+1):
        for j in range(1,h_resi.GetNbinsY()+1):
            val = h_resi.GetBinContent(i,j)
            h_resi.SetBinContent(i,j,val/h.GetBinError(i,j))
    h_resi = set2DHisto(h_resi,xTitle,yTitle,zTitle)
    absMax = max(abs(h_resi.GetMinimum()),abs(h_resi.GetMaximum()))
    setFFColors(h_resi,-1.5*absMax,1.5*absMax)
    h_resi.Draw("colztext")
    xLines, yLines = getGrayLines(x,y)
    
    [xLine.Draw("l") for xLine in xLines]
    [yLine.Draw("l") for yLine in yLines]

    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.91,"CMS simulation")
    l.DrawLatex(0.65,0.91,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.DrawLatex(0.2,0.85,boxLabel)
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
        

def get3DHistoFrom1D(h1D,x,y,z,name):   
    h3D = rt.TH3D(name,name,len(x)-1,x,len(y)-1,y,len(z)-1,z)

    iBinX=-1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                h3D.SetBinContent(i,j,k,h1D.GetBinContent(iBinX+1))
                h3D.SetBinError(i,j,k,h1D.GetBinError(iBinX+1))
    return h3D

def Gamma(a, x):
    return rt.TMath.Gamma(a) * rt.Math.inc_gamma_c(a,x)

def Gfun(x, y, X0, Y0, B, N):
    return Gamma(N,B*N*rt.TMath.Power((x-X0)*(y-Y0),1/N))

def getBinEvents(i, j, k, x, y, z, workspace):
    errorFlag = False
    bkg = "TTj%ib"%z[k-1]
    if z[k-1]==3:
        bkgShape = "TTj2b"
    else:
        bkgShape = "TTj%ib"%z[k-1]
        
    B = workspace.var("b_%s_%s"%(bkgShape,box)).getVal()
    N = workspace.var("n_%s_%s"%(bkgShape,box)).getVal()
    X0 = workspace.var("MR0_%s_%s"%(bkgShape,box)).getVal()
    Y0 = workspace.var("R0_%s_%s"%(bkgShape,box)).getVal()
    NTOT = workspace.var("Ntot_%s_%s"%(bkg,box)).getVal()

    xmin  = x[0]
    xmax  = x[-1]
    ymin  = y[0]
    ymax  = y[-1]
    total_integral = (N/rt.TMath.Power(B*N,N))*(Gfun(xmin,ymin,X0,Y0,B,N)-Gfun(xmin,ymax,X0,Y0,B,N)-Gfun(xmax,ymin,X0,Y0,B,N)+Gfun(xmax,ymax,X0,Y0,B,N))

    xmin  = x[i-1]
    xmax  = x[i]
    ymin  = y[j-1]
    ymax  = y[j]
    integral = (N/rt.TMath.Power(B*N,N))*(Gfun(xmin,ymin,X0,Y0,B,N)-Gfun(xmin,ymax,X0,Y0,B,N)-Gfun(xmax,ymin,X0,Y0,B,N)+Gfun(xmax,ymax,X0,Y0,B,N))

    bin_events =  NTOT*integral/total_integral

    if bin_events <= 0:
        errorFlag = True
        print "\nERROR: bin razor pdf integral =", integral
        print "\nERROR: total razor pdf integral =", total_integral
        return 0., errorFlag
    return bin_events, errorFlag



if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store results")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-i','--inputFitFile',dest="inputFitFile", default="BinnedFitResults_MultiJet.root",type="string",
                  help="intput fit file")
    parser.add_option('-t','--inputToyFile',dest="inputToyFile", default="toys_Bayes_MultiJet.root",type="string",
                  help="intput fit file")
    
    (options,args) = parser.parse_args()
     
    box = options.box
    lumi = options.lumi
    cfg = Config.Config(options.config)

    inputFitFile = rt.TFile.Open(options.inputFitFile,"read")

    w = inputFitFile.Get("w"+box)

    w.Print("v")
        
    th1x = w.var('th1x')
    dataHist = w.data("data_obs")

    setStyle()
    
    extRazorPdf = w.pdf('extRazorPdf')
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)

    
    xFine = array('d', [x[0]+i*(x[-1]-x[0])/100. for i in range(0,101)]) # MR binning fine
    yFine = array('d', [y[0]+i*(y[-1]-y[0])/100. for i in range(0,101)]) # Rsq binning fine
    zFine = array('d', cfg.getBinning(box)[2]) # nBtag binning fine
    nBinsFine = (len(xFine)-1)*(len(yFine)-1)*(len(zFine)-1)
    
    
    th1x.setBins(nBins)

    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())
    
    
    rt.TH1D.SetDefaultSumw2()
    rt.TH2D.SetDefaultSumw2()
    rt.TH3D.SetDefaultSumw2()

    h_th1x = asimov.createHistogram('h_th1x',th1x)
    h_data_th1x = dataHist.createHistogram('h_data_th1x',th1x)
    
    h_data_nBtagRsqMR_fine = rt.TH3D("h_data_nBtagRsqMR_fine","h_data_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    h_nBtagRsqMR_fine = rt.TH3D("h_nBtagRsqMR_fine","h_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    w.data("RMRTree").fillHistogram(h_data_nBtagRsqMR_fine,rt.RooArgList(w.var("MR"),w.var("Rsq"),w.var("nBtag")))
    for i in range(1,len(xFine)):
        for j in range(1,len(yFine)):
            for k in range(1,len(zFine)):
                value, errorFlag= getBinEvents(i,j,k,xFine,yFine,zFine,w)
                if not errorFlag:
                    h_nBtagRsqMR_fine.SetBinContent(i,j,k,value)

                    
    h_data_RsqMR_fine = h_data_nBtagRsqMR_fine.Project3D("yxe")
    h_RsqMR_fine = h_nBtagRsqMR_fine.Project3D("yxe")
    
    h_data_nBtagRsqMR = get3DHistoFrom1D(h_data_th1x,x,y,z,"h_data_nBtagRsqMR")
    h_nBtagRsqMR = get3DHistoFrom1D(h_th1x,x,y,z,"h_nBtagRsqMR")
    
        
    h_data_RsqMR = h_data_nBtagRsqMR.Project3D("yxe")
    h_data_MR = h_data_nBtagRsqMR.Project3D("xe")
    h_data_Rsq = h_data_nBtagRsqMR.Project3D("ye")
    h_RsqMR = h_nBtagRsqMR.Project3D("yxe")
    h_MR = h_nBtagRsqMR.Project3D("xe")
    h_Rsq = h_nBtagRsqMR.Project3D("ye")


    h_MR_slices = []
    h_MR_integrals = []
    h_data_MR_slices = []
    h_data_MR_integrals = []
    h_MR_slice_labels = []
    h_MR_integral_labels = []
    
    h_MR_slice_components = []
    h_data_MR_slice_components = []
    h_MR_slice_component_labels = []
    for j in range(1,len(y)):
        h_MR_slices.append(h_nBtagRsqMR.ProjectionX("MR_%.2fRsq%.2f"%(y[j-1],y[j]),j,j,0,-1,""))
        h_MR_integrals.append(h_nBtagRsqMR.ProjectionX("MR_Rsq%.2f"%(y[j-1]),j,len(y)-1,0,-1,""))
        h_data_MR_slices.append(h_data_nBtagRsqMR.ProjectionX("MR_data_%.2fRsq%.2f"%(y[j-1],y[j]),j,j,0,-1,""))
        h_data_MR_integrals.append(h_data_nBtagRsqMR.ProjectionX("MR_data_Rsq%.2f"%(y[j-1]),j,len(y)-1,0,-1,""))
        h_MR_slice_labels.append("%.2f #leq R^{2} < %.2f"%(y[j-1],y[j]))
        h_MR_integral_labels.append("R^{2} #geq %.2f"%(y[j-1]))
        if len(z)>1:
            h_MR_comp = []
            h_data_MR_comp = []
            h_label_comp = []
            for k in range(1,len(z)):
                h_MR_comp.append(h_nBtagRsqMR.ProjectionX("MR_%ibtag_%.2fRsq%.2f"%(z[k-1],y[j-1],y[j]),j,j,k,k,""))
                h_data_MR_comp.append(h_data_nBtagRsqMR.ProjectionX("MR_data_%ibtag_%.2fRsq%.2f"%(z[k-1],y[j-1],y[j]),j,j,k,k,""))                
                if z[k-1]==3 and z[-1]==4:
                    h_label_comp.append("#geq %i b-tag, %.2f #leq R^{2} < %.2f " % (z[k-1],y[j-1],y[j]) )
                if z[k-1]==1 and z[-1]==2 and box in ['MuEle','MuMu','EleEle']:                
                    h_label_comp.append("#geq %i b-tag, %.2f #leq R^{2} < %.2f " % (z[k-1],y[j-1],y[j]) )
                else:            
                    h_label_comp.append("%i b-tag, %.2f #leq R^{2} < %.2f " % (z[k-1],y[j-1],y[j]) )
            h_MR_slice_components.append(h_MR_comp)
            h_data_MR_slice_components.append(h_data_MR_comp)
            h_MR_slice_component_labels.append(h_label_comp)
            

        
    h_Rsq_slices = []
    h_Rsq_integrals = []
    h_data_Rsq_slices = []
    h_data_Rsq_integrals = []
    h_Rsq_slice_labels = []
    h_Rsq_integral_labels = []
    
    h_Rsq_slice_components = []
    h_data_Rsq_slice_components = []
    h_Rsq_slice_component_labels = []
    for i in range(1,len(x)):
        h_Rsq_slices.append(h_nBtagRsqMR.ProjectionY("Rsq_%iMR%i"%(x[i-1],x[i]),i,i,0,-1,""))
        h_Rsq_integrals.append(h_nBtagRsqMR.ProjectionY("Rsq_MR%i"%(x[i-1]),i,len(x)-1,0,-1,""))
        h_data_Rsq_slices.append(h_data_nBtagRsqMR.ProjectionY("Rsq_data_%iMR%i"%(x[i-1],x[i]),i,i,0,-1,""))
        h_data_Rsq_integrals.append(h_data_nBtagRsqMR.ProjectionY("Rsq_data_MR%i"%(x[i-1]),i,len(x)-1,0,-1,""))
        h_Rsq_slice_labels.append("%i #leq M_{R} < %i"%(x[i-1],x[i]))
        h_Rsq_integral_labels.append("M_{R} #geq %i"%(x[i-1]))
        if len(z)>1:
            h_Rsq_comp = []
            h_data_Rsq_comp = []
            h_label_comp = []
            for k in range(1,len(z)):
                h_Rsq_comp.append(h_nBtagRsqMR.ProjectionY("Rsq_%ibtag_%iMR%i"%(z[k-1],x[i-1],x[i]),i,i,k,k,""))
                h_data_Rsq_comp.append(h_data_nBtagRsqMR.ProjectionY("Rsq_data_%ibtag_%iMR%i"%(z[k-1],x[i-1],x[i]),i,i,k,k,""))                
                if z[k-1]==3 and z[-1]==4:
                    h_label_comp.append("#geq %i b-tag, %i #leq M_{R} < %i " % (z[k-1],x[i-1],x[i]) )
                if z[k-1]==1 and z[-1]==2 and box in ['MuEle','MuMu','EleEle']:                
                    h_label_comp.append("#geq %i b-tag, %i #leq M_{R} < %i " % (z[k-1],x[i-1],x[i]) )
                else:            
                    h_label_comp.append("%i b-tag, %i #leq M_{R} < %i " % (z[k-1],x[i-1],x[i]) )
            h_Rsq_slice_components.append(h_Rsq_comp)
            h_data_Rsq_slice_components.append(h_data_Rsq_comp)
            h_Rsq_slice_component_labels.append(h_label_comp)

                        
    h_MR_components = []
    h_Rsq_components = []
    h_data_MR_components = []
    h_data_Rsq_components = []
    h_labels = []        
    h_colors = [rt.kOrange,rt.kViolet,rt.kRed,rt.kGreen]
    if len(z)>1:
        for k in range(1,len(z)):
            h_MR_components.append(h_nBtagRsqMR.ProjectionX("MR_%ibtag"%z[k-1],0,-1,k,k,""))
            h_Rsq_components.append(h_nBtagRsqMR.ProjectionY("Rsq_%ibtag"%z[k-1],0,-1,k,k,""))
            h_data_MR_components.append(h_data_nBtagRsqMR.ProjectionX("MR_data_%ibtag"%z[k-1],0,-1,k,k,""))
            h_data_Rsq_components.append(h_data_nBtagRsqMR.ProjectionY("Rsq_data_%ibtag"%z[k-1],0,-1,k,k,""))
            if z[k-1]==3 and z[-1]==4:
                h_labels.append("#geq %i b-tag" % z[k-1] )
            if z[k-1]==1 and z[-1]==2 and box in ['MuEle','MuMu','EleEle']:                
                h_labels.append("#geq %i b-tag" % z[k-1] )
            else:            
                h_labels.append("%i b-tag" % z[k-1] )


    btagLabel = ""
    if z[-1] == z[0]+1 and z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    elif z[-1] == z[0]+1:
        btagLabel = "%i b-tag" % z[0]
    elif z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    elif z[-1]==2 and box in ['MuEle','MuMu','EleEle']:
        btagLabel = "#geq %i b-tag" % z[0]        
    else:
        btagLabel = "%i-%i b-tag" % (z[0],z[-2])

    lumiLabel = "%i fb^{-1} (13 TeV)" % (lumi/1000.)
    boxLabel = "razor %s %s" % (box,btagLabel)

    
    c = rt.TCanvas('c','c',500,400)
    print1DProj(c,h_MR,h_data_MR,options.outDir+"/h_MR_%s.pdf"%box,"M_{R} [GeV]","Events",lumiLabel,boxLabel,None,h_MR_components,h_colors,h_labels)
    print1DProj(c,h_Rsq,h_data_Rsq,options.outDir+"/h_Rsq_%s.pdf"%box,"R^{2}","Events",lumiLabel,boxLabel,None,h_Rsq_components,h_colors,h_labels)

    
    newBoxLabel = "razor %s %s" % (box,btagLabel)
    more_colors = [rt.kBlack,rt.kBlue]
    more_colors.extend(h_colors)
    print1DSlice(c,h_MR_integrals,h_data_MR_integrals,options.outDir+"/h_MR_slicesRsq_%s.pdf"%(box),"M_{R} [GeV]","Events",lumiLabel,newBoxLabel,None,more_colors,h_MR_integral_labels)
    
    newBoxLabel = "razor %s %s" % (box,btagLabel)
    more_colors = [rt.kBlack,rt.kBlue]
    more_colors.extend(h_colors)
    more_colors.extend([rt.kMagenta,rt.kGray])
    print1DSlice(c,h_Rsq_integrals,h_data_Rsq_integrals,options.outDir+"/h_Rsq_slicesMR_%s.pdf"%(box),"R^{2}","Events",lumiLabel,newBoxLabel,None,more_colors,h_Rsq_integral_labels)

    print2DResiduals(c,h_RsqMR,h_data_RsqMR,options.outDir+"/h_RsqMR_residuals_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Residuals (Sim. Data - Fit)",lumiLabel,boxLabel,x,y)
    print2DPercentDiff(c,h_RsqMR,h_data_RsqMR,options.outDir+"/h_RsqMR_percentdiff_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Percent Diff. (Sim. Data - Fit)/Fit",lumiLabel,boxLabel,x,y)
    print2DNSigma(c,h_RsqMR,h_data_RsqMR,options.outDir+"/h_RsqMR_nsigma_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "n#sigma (Sim. Data - Fit)/sqrt(Fit)",lumiLabel,boxLabel,x,y)
    print2DScatter(c,h_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Fit",lumiLabel,boxLabel,x,y,h_data_RsqMR_fine.GetMinimum(),h_data_RsqMR_fine.GetMaximum())
    print2DScatter(c,h_data_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_data_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Sim. Data",lumiLabel,boxLabel,x,y,h_data_RsqMR_fine.GetMinimum(),h_data_RsqMR_fine.GetMaximum())

    #sys.exit()
    
    for k in range(0,len(z)-1):
        newBoxLabel = "razor %s %s"%(box,h_labels[k])
        print1DProj(c,h_MR_components[k],h_data_MR_components[k],options.outDir+"/h_MR_%ibtag_%s.pdf"%(z[k],box),"M_{R} [GeV]","Events",lumiLabel,newBoxLabel)
        print1DProj(c,h_Rsq_components[k],h_data_Rsq_components[k],options.outDir+"/h_Rsq_%ibtag_%s.pdf"%(z[k],box),"R^{2}","Events",lumiLabel,newBoxLabel)

        
    for j in range(0,len(y)-1):
        newBoxLabel = "razor %s %s"%(box,h_MR_slice_labels[j])
        print1DProj(c,h_MR_slices[j],h_data_MR_slices[j],options.outDir+"/h_MR_%.2fRsq%.2f_%s.pdf"%(y[j],y[j+1],box),"M_{R} [GeV]","Events",lumiLabel,newBoxLabel,None,h_MR_slice_components[j],h_colors,h_labels)
        for k in range(0,len(z)-1):
            newBoxLabel = "razor %s %s"%(box,h_MR_slice_component_labels[j][k])
            print1DProj(c,h_MR_slice_components[j][k],h_data_MR_slice_components[j][k],options.outDir+"/h_MR_%ibtag_%.2fRsq%.2f_%s.pdf"%(z[k],y[j],y[j+1],box),"M_{R} [GeV]","Events",lumiLabel,newBoxLabel)
        
    for i in range(0,len(x)-1):
        newBoxLabel = "razor %s %s"%(box,h_Rsq_slice_labels[i])
        print1DProj(c,h_Rsq_slices[i],h_data_Rsq_slices[i],options.outDir+"/h_Rsq_%iMR%i_%s.pdf"%(x[i],x[i+1],box),"R^{2}","Events",lumiLabel,newBoxLabel,None,h_Rsq_slice_components[i],h_colors,h_labels)
        for k in range(0,len(z)-1):
            newBoxLabel = "razor %s %s"%(box,h_Rsq_slice_component_labels[i][k])
            print1DProj(c,h_Rsq_slice_components[i][k],h_data_Rsq_slice_components[i][k],options.outDir+"/h_Rsq_%ibtag_%iMR%i_%s.pdf"%(z[k],x[i],x[i+1],box),"R^{2}","Events",lumiLabel,newBoxLabel)

            
