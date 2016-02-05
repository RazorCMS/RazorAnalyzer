from optparse import OptionParser
import ROOT as rt
import sys
import rootTools
import glob
from math import *
import os
from array import *
import numpy as np
from scipy.interpolate import Rbf
import itertools
from GChiPairs import gchipairs

toFix = []
def interpolate2D(hist,epsilon=1,smooth=0,diagonalOffset=0):
    x = array('d',[])
    y = array('d',[])
    z = array('d',[])
    
    binWidth = float(hist.GetXaxis().GetBinWidth(1))
    
    for i in range(1, hist.GetNbinsX()+1):
        for j in range(1, hist.GetNbinsY()+1):
            if hist.GetBinContent(i,j)>0.:
                x.append(hist.GetXaxis().GetBinCenter(i))
                y.append(hist.GetYaxis().GetBinCenter(j))
                z.append(rt.TMath.Log(hist.GetBinContent(i,j)))
                #z.append(hist.GetBinContent(i,j))

    mgMin = hist.GetXaxis().GetBinCenter(1)
    mgMax = hist.GetXaxis().GetBinCenter(hist.GetNbinsX())#+hist.GetXaxis().GetBinWidth(hist.GetNbinsX())
    mchiMin = hist.GetYaxis().GetBinCenter(1)
    mchiMax = hist.GetYaxis().GetBinCenter(hist.GetNbinsY())#+hist.GetYaxis().GetBinWidth(hist.GetNbinsY())
    
    myX = np.linspace(mgMin, mgMax,int((mgMax-mgMin)/binWidth+1))
    myY = np.linspace(mchiMin, mchiMax, int((mchiMax-mchiMin)/binWidth+1))
    myXI, myYI = np.meshgrid(myX,myY)

    rbf = Rbf(x, y, z,function='multiquadric', epsilon=epsilon,smooth=smooth)
    myZI = rbf(myXI, myYI)

    for i in range(1, hist.GetNbinsX()+1):
        for j in range(1, hist.GetNbinsY()+1):
            xLow = hist.GetXaxis().GetBinCenter(i)
            yLow = hist.GetYaxis().GetBinCenter(j)
            if xLow >= yLow+diagonalOffset-binWidth/2.:
                hist.SetBinContent(i,j,rt.TMath.Exp(myZI[j-1][i-1]))
    return hist


def fix_hist_byhand(hist, model, box, clsType):
    if 'Obs' in clsType:
        for (mg,mchi) in gchipairs(model):            
            #obs = xsecUL['Obs'].GetBinContent(xsecUL['Obs'].FindBin(mg,mchi))
            obs = hist.GetBinContent(hist.FindBin(mg,mchi))
            exp = xsecUL['Exp'].GetBinContent(xsecUL['Exp'].FindBin(mg,mchi))
            expPlus2 = xsecUL['ExpPlus2'].GetBinContent(xsecUL['ExpPlus2'].FindBin(mg,mchi))
            expPlus = xsecUL['ExpPlus'].GetBinContent(xsecUL['ExpPlus'].FindBin(mg,mchi))
            expMinus = xsecUL['ExpMinus'].GetBinContent(xsecUL['ExpMinus'].FindBin(mg,mchi))
            expMinus2 = xsecUL['ExpMinus2'].GetBinContent(xsecUL['ExpMinus2'].FindBin(mg,mchi))
            if hist.GetBinContent(hist.FindBin(mg,mchi))==0:
                if (mg,mchi) not in toFix:                
                    #print "empty:", mg, mchi
                    toFix.append((mg,mchi))
            elif obs<expPlus2 or obs>expMinus2:
                hist.SetBinContent(hist.FindBin(mg,mchi),0)
                if (mg,mchi) not in toFix:
                    #print "outside 1s:", mg, mchi
                    toFix.append((mg,mchi))

        if model == "T1tttt" and box == "MultiJet":
             hist.SetBinContent(hist.FindBin(975,700),0)            
        ## if model == "T1tttt" and box == "MuMultiJet_EleMultiJet":
        ##     hist.SetBinContent(hist.FindBin(950,600),0)
        ##     hist.SetBinContent(hist.FindBin(1400,700),0)
        ##     hist.SetBinContent(hist.FindBin(1500,300),0)
        ##     hist.SetBinContent(hist.FindBin(1800,300),0)
        ##     hist.SetBinContent(hist.FindBin(1850,50),0)
        ##     hist.SetBinContent(hist.FindBin(850,300),0)
        ##     hist.SetBinContent(hist.FindBin(1900,850),0)
        ## if model == "T1tttt" and box == "MuMultiJet_EleMultiJet_MultiJet":
        ##     hist.SetBinContent(hist.FindBin(1700,0),0)
        ##     hist.SetBinContent(hist.FindBin(700,200),0)
        ##     hist.SetBinContent(hist.FindBin(900,200),0)
        ##     hist.SetBinContent(hist.FindBin(1100,400),0)
        ##     hist.SetBinContent(hist.FindBin(1200,100),0)
        ##     hist.SetBinContent(hist.FindBin(1200,400),0)
        ##     hist.SetBinContent(hist.FindBin(1200,750),0)
        ##     hist.SetBinContent(hist.FindBin(1250,875),0)
                
                        
def set_palette(name="default", ncontours=255):
    # For the canvas:
    rt.gStyle.SetCanvasBorderMode(0)
    rt.gStyle.SetCanvasColor(rt.kWhite)
    rt.gStyle.SetCanvasDefH(400) #Height of canvas
    rt.gStyle.SetCanvasDefW(500) #Width of canvas
    rt.gStyle.SetCanvasDefX(0)   #POsition on screen
    rt.gStyle.SetCanvasDefY(0)
	
    # For the Pad:
    rt.gStyle.SetPadBorderMode(0)
    # rt.gStyle.SetPadBorderSize(Width_t size = 1)
    rt.gStyle.SetPadColor(rt.kWhite)
    rt.gStyle.SetPadGridX(False)
    rt.gStyle.SetPadGridY(False)
    rt.gStyle.SetGridColor(0)
    rt.gStyle.SetGridStyle(3)
    rt.gStyle.SetGridWidth(1)
	
    # For the frame:
    rt.gStyle.SetFrameBorderMode(0)
    rt.gStyle.SetFrameBorderSize(1)
    rt.gStyle.SetFrameFillColor(0)
    rt.gStyle.SetFrameFillStyle(0)
    rt.gStyle.SetFrameLineColor(1)
    rt.gStyle.SetFrameLineStyle(1)
    rt.gStyle.SetFrameLineWidth(1)
    
    # set the paper & margin sizes
    rt.gStyle.SetPaperSize(20,26)
    rt.gStyle.SetPadTopMargin(0.085)
    rt.gStyle.SetPadRightMargin(0.15)
    rt.gStyle.SetPadBottomMargin(0.15)
    rt.gStyle.SetPadLeftMargin(0.17)
    
    # use large Times-Roman fonts
    rt.gStyle.SetTitleFont(42,"xyz")  # set the all 3 axes title font
    rt.gStyle.SetTitleFont(42," ")    # set the pad title font
    rt.gStyle.SetTitleSize(0.06,"xyz") # set the 3 axes title size
    rt.gStyle.SetTitleSize(0.06," ")   # set the pad title size
    rt.gStyle.SetLabelFont(42,"xyz")
    rt.gStyle.SetLabelSize(0.05,"xyz")
    rt.gStyle.SetLabelColor(1,"xyz")
    rt.gStyle.SetTextFont(42)
    rt.gStyle.SetTextSize(0.08)
    rt.gStyle.SetStatFont(42)
    
    # use bold lines and markers
    rt.gStyle.SetMarkerStyle(8)
    #rt.gStyle.SetHistLineWidth(1.85)
    rt.gStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    
    #..Get rid of X error bars
    rt.gStyle.SetErrorX(0.001)
    
    # do not display any of the standard histogram decorations
    rt.gStyle.SetOptTitle(1)
    rt.gStyle.SetOptStat(0)
    #rt.gStyle.SetOptFit(11111111)
    rt.gStyle.SetOptFit(0)
    
    # put tick marks on top and RHS of plots
    rt.gStyle.SetPadTickX(1)
    rt.gStyle.SetPadTickY(1)
    if name == "gray" or name == "grayscale":
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [1.00, 0.95, 0.95, 0.65, 0.15]
        green = [1.00, 0.85, 0.7, 0.5, 0.3]
        blue  = [0.95, 0.6, 0.3, 0.45, 0.65]
    elif name == "chris":
	stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
	red =   [ 1.0,   0.95,  0.95,  0.65,   0.15 ]
	green = [ 1.0,  0.85, 0.7, 0.5,  0.3 ]
	blue =  [ 0.95, 0.6 , 0.3,  0.45, 0.65 ]
    elif name == "blue":
	stops = [ 0.0, 0.5, 1.0]
	red =   [ 0.0, 1.0, 0.0]
	green = [ 0.0, 0.0, 0.0]
	blue =  [ 1.0, 0.0, 0.0]
    else:
        # default palette, looks cool
        stops = [0.00, 0.34, 0.61, 0.84, 1.00]
        red   = [0.00, 0.00, 0.87, 1.00, 0.51]
        green = [0.00, 0.81, 1.00, 0.20, 0.00]
        blue  = [0.51, 1.00, 0.12, 0.00, 0.00]
        #from PlotsSMS/python/smsPlotXSEC.py
        stops = [0.00, 0.20, 0.70, 0.90, 1.00]
        red   = [0.00, 0.00, 1.00, 1.00, 1.00]
        green = [0.00, 1.00, 1.00, 0.30, 0.00]
        blue  = [1.00, 1.00, 0.00, 0.20, 0.00]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    
    npoints = len(s)
    rt.TColor.CreateGradientColorTable(npoints, s, r, g, b, ncontours)
    rt.gStyle.SetNumberContours(ncontours)

def getLogHist(hist):
    logHist = hist.Clone("log"+hist.GetName())
    for i in xrange(1,hist.GetNbinsX()+1):
        for j in xrange(1,hist.GetNbinsY()+1):
            if hist.GetBinContent(i,j) > 0.: logHist.SetBinContent(i,j,rt.TMath.Log(hist.GetBinContent(i,j)))
    return logHist

def getExpHist(logHist):
    hist = logHist.Clone("exp"+logHist.GetName())
    for i in xrange(1,hist.GetNbinsX()+1):
        for j in xrange(1,hist.GetNbinsY()+1):
            if logHist.GetBinContent(i,j) != 0.: hist.SetBinContent(i,j,rt.TMath.Exp(logHist.GetBinContent(i,j)))
    return hist

def getModelSettings(model):
    if model=="T1bbbb":
        mgMin = 600.-12.5
        mgMax = 2000.+12.5
        mchiMin = 0.-12.5
        mchiMax = 1450.+12.5 
        binWidth = 25
        nRebins = 0
        xsecMin = 1.e-3
        xsecMax = 10.
        diagonalOffset = 25+12.5
        smoothing = 50
    elif model=="T1tttt":
        mgMin = 600.-12.5
        mgMax = 2000.+12.5
        mchiMin = 0.-12.5
        mchiMax = 1450.+12.5 
        binWidth = 25
        nRebins = 0
        xsecMin = 1.e-3
        xsecMax = 10.
        diagonalOffset = 225+12.5
        smoothing = 200
    elif model=="T1qqqq":
        mgMin = 600.-12.5
        mgMax = 2000.+12.5
        mchiMin = 0.-12.5
        mchiMax = 1250.+12.5
        binWidth = 25
        nRebins = 0
        xsecMin = 1.e-3
        xsecMax = 10.
        diagonalOffset = 25+12.5
        smoothing = 50
    return mgMin, mgMax, mchiMin, mchiMax, binWidth, nRebins, xsecMin, xsecMax, diagonalOffset, smoothing

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")
    parser.add_option('--toys',dest="doHybridNew",default=False,action='store_true',
                  help="for toys instead of asymptotic")
    parser.add_option('--xsec-file',dest="refXsecFile",default="./data/gluino13TeV.txt",type="string",
                  help="Input directory")
    
    (options,args) = parser.parse_args()
    
    box = options.box
    model = options.model
    directory = options.outDir

    
    refXsecFile = options.refXsecFile
    doHybridNew = options.doHybridNew
                
    set_palette("rainbow",255)
    rt.gStyle.SetOptStat(0)
    rt.gROOT.ProcessLine(".L macros/swissCrossInterpolate.h++")
    rt.gSystem.Load("macros/swissCrossInterpolate_h.so")

    mgMin, mgMax, mchiMin, mchiMax, binWidth, nRebins, xsecMin, xsecMax, diagonalOffset, smoothing = getModelSettings(model)

    if doHybridNew:
        xsecFile = rt.TFile.Open("%s/xsecUL_HybridNew_%s.root"%(directory,box))
    else: 
        xsecFile = rt.TFile.Open("%s/xsecUL_Asymptotic_%s.root"%(directory,box))
    xsecTree = xsecFile.Get("xsecTree")
    xsecGluino =  rt.TH2D("xsecGluino","xsecGluino",int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
    xsecGluinoPlus =  rt.TH2D("xsecGluinoPlus","xsecGluinoPlus",int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
    xsecGluinoMinus =  rt.TH2D("xsecGluinoMinus","xsecGluinoMinus",int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
    
    clsTypes = ["Exp","ExpMinus","ExpMinus2","ExpPlus","ExpPlus2","Obs","ObsMinus","ObsPlus"]

    smooth = {}
    for clsType in clsTypes:
        smooth[clsType] = smoothing
    
    titleMap = {"Exp":"Expected","ExpMinus":"Expected-1#sigma","ExpPlus":"Expected+1#sigma",
                "ExpMinus2":"Expected-2#sigma","ExpPlus2":"Expected+2#sigma",
                "ObsMinus":"Observed-1#sigma", "ObsPlus":"Observed+1#sigma","Obs":"Observed"}
    whichCLsVar = {"Obs":"xsecULObs_%s"%(box),"ObsPlus":"xsecULObs_%s"%(box),"ObsMinus":"xsecULObs_%s"%(box),
                   "Exp":"xsecULExp_%s"%(box),"ExpPlus":"xsecULExpMinus_%s"%(box),"ExpMinus":"xsecULExpPlus_%s"%(box),
                   "ExpPlus2":"xsecULExpMinus2_%s"%(box),"ExpMinus2":"xsecULExpPlus2_%s"%(box)}
                   
    xsecUL = {}
    logXsecUL = {}
    rebinXsecUL = {}
    subXsecUL = {}

    contourFinal = {}

        
    
    subboxes = box.split("_")
    
    for clsType in clsTypes:
        xsecUL[clsType] = rt.TH2D("xsecUL_%s"%clsType,"xsecUL_%s"%clsType,int((mgMax-mgMin)/binWidth),mgMin, mgMax,int((mchiMax-mchiMin)/binWidth), mchiMin, mchiMax)
        xsecTree.Project("xsecUL_%s"%clsType,"mchi:mg",whichCLsVar[clsType])


        fix_hist_byhand(xsecUL[clsType],model,box,clsType)
        
            
        print "INFO: doing interpolation for %s"%(clsType)
        
        # do swiss cross average in real domain
        rebinXsecUL[clsType] = rt.swissCrossInterpolate(xsecUL[clsType],"NE")
        #rebinXsecUL[clsType] = xsecUL[clsType].Clone()

        # do scipy multi-quadratic interpolation in log domain
        rebinXsecUL[clsType] = interpolate2D(rebinXsecUL[clsType],epsilon=5,smooth=smooth[clsType],diagonalOffset=diagonalOffset)

        # do swiss cross rebin + average in real domain (should be log??)
        for i in xrange(0,nRebins):
            rebinXsecUL[clsType] = rt.swissCrossRebin(rebinXsecUL[clsType],"NE")

        # only for display purposes of underlying heat map: do swiss cross average then scipy interpolation 
        xsecUL[clsType] = rt.swissCrossInterpolate(xsecUL[clsType],"NE")
        xsecUL[clsType] = interpolate2D(xsecUL[clsType], epsilon=5,smooth=smooth[clsType],diagonalOffset=diagonalOffset)

        # fix axes
        xsecUL[clsType].GetXaxis().SetRangeUser(xsecUL[clsType].GetXaxis().GetBinCenter(1),xsecUL[clsType].GetXaxis().GetBinCenter(xsecUL[clsType].GetNbinsX()))
        xsecUL[clsType].GetYaxis().SetRangeUser(xsecUL[clsType].GetYaxis().GetBinCenter(1),xsecUL[clsType].GetYaxis().GetBinCenter(xsecUL[clsType].GetNbinsY()))
         
        
    thyXsec = {}
    thyXsecErr = {}
    if refXsecFile is not None:
        print "INFO: Input ref xsec file!"
        for mg in range(600,2025,25):
            for line in open(refXsecFile,'r'):
                line = line.replace('\n','')
                if str(mg)==line.split(',')[0]:
                    thyXsec[mg] = float(line.split(',')[1]) #pb
                    thyXsecErr[mg] = 0.01*float(line.split(',')[2])
    else: 
        print "ERROR: no xsec file; exiting"
        sys.exit()  
    

    # now rebin xsecGluino the correct number of times
    for i in xrange(0,nRebins):
        xsecGluino = rt.swissCrossRebin(xsecGluino,"NE")
        xsecGluinoPlus = rt.swissCrossRebin(xsecGluinoPlus,"NE")
        xsecGluinoMinus = rt.swissCrossRebin(xsecGluinoMinus,"NE")
    
    
    for i in xrange(1,xsecGluino.GetNbinsX()+1):
        xLow = xsecGluino.GetXaxis().GetBinCenter(i)
        xsecVal = thyXsec[int(xLow)]
        xsecErr =  thyXsecErr[int(xLow)]
        for j in xrange(1,xsecGluino.GetNbinsY()+1):
            yLow = xsecGluino.GetYaxis().GetBinCenter(j)
            if xLow >= yLow+diagonalOffset and xLow <= mgMax-binWidth/2:
                xsecGluino.SetBinContent(i,j,xsecVal)
                xsecGluinoPlus.SetBinContent(i,j,xsecVal*(1+xsecErr))
                xsecGluinoMinus.SetBinContent(i,j,xsecVal*(1-xsecErr))
                
    c = rt.TCanvas("c","c",500,500)
    
    for clsType in clsTypes:
        xsecUL[clsType].SetName("xsecUL_%s_%s_%s"%(clsType,model,box))
        xsecUL[clsType].SetTitle("%s %s, %s #sigma #times Branching Fraction"%(model,box.replace("_","+"),titleMap[clsType]))
        xsecUL[clsType].SetTitleOffset(1.5)
        xsecUL[clsType].SetXTitle("gluino Mass [GeV]")
        xsecUL[clsType].GetXaxis().SetTitleOffset(0.95)
        xsecUL[clsType].SetYTitle("LSP Mass [GeV]")
        xsecUL[clsType].GetYaxis().SetTitleOffset(1.25)

        # subtract the predicted xsec
        subXsecUL[clsType] = rebinXsecUL[clsType].Clone()
        if clsType=="ObsMinus":
            subXsecUL[clsType].Add(xsecGluinoMinus,-1)
        elif clsType=="ObsPlus":
            subXsecUL[clsType].Add(xsecGluinoPlus,-1)
        else:
            subXsecUL[clsType].Add(xsecGluino,-1)

        contours = array('d',[0.0])
        subXsecUL[clsType].SetContour(1,contours)
        
        xsecUL[clsType].SetMinimum(xsecMin)
        xsecUL[clsType].SetMaximum(xsecMax)
        rebinXsecUL[clsType].SetMinimum(xsecMin)
        rebinXsecUL[clsType].SetMaximum(xsecMax)
        subXsecUL[clsType].SetMaximum(1.)
        subXsecUL[clsType].SetMinimum(-1.)
        
        c.SetLogz(0)
        subXsecUL[clsType].Draw("CONT Z LIST")
        c.Update()
        
        conts = rt.gROOT.GetListOfSpecials().FindObject("contours")

        xsecUL[clsType].Draw("COLZ")
        #subXsecUL[clsType].Draw("COLZ")
        
        contour0 = conts.At(0)
        curv = contour0.First()
        finalcurv = rt.TGraph(1)
        try:
            curv.SetLineWidth(3)
            curv.SetLineColor(rt.kBlack)
            curv.Draw("lsame")
            finalcurv = curv.Clone()
            maxN = curv.GetN()
        except AttributeError:
            print "ERROR: no curve drawn for box=%s, clsType=%s -- no limit "%(box, clsType)
        
        for i in xrange(1, contour0.GetSize()):
            curv = contour0.After(curv)
            curv.SetLineWidth(3)
            curv.SetLineColor(rt.kBlack)
            curv.Draw("lsame")
            if curv.GetN()>maxN:
                maxN = curv.GetN()
                finalcurv = curv.Clone()

        contourFinal[clsType] = finalcurv

        contourFinal[clsType].SetName("%s_%s_%s"%(clsType,model,box))
        
        c.SetLogz(1)
        c.Print("%s/%s_INTERP_%s_%s.pdf"%(directory,model,box,clsType))

    outFile = rt.TFile.Open("%s/%s_%s_results.root"%(directory,model,box),"recreate")
    for clsType in clsTypes:
        contourFinal[clsType].Write()
        xsecUL[clsType].Write()

    smoothOutFile = rt.TFile.Open("%s/smoothXsecUL_%s_%s.root"%(directory,model,box), "recreate")
    
    smoothXsecTree = rt.TTree("smoothXsecTree", "smoothXsecTree")
    myStructCmd = "struct MyStruct{Double_t mg;Double_t mchi;"
    ixsecUL = 0
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+0)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+1)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+2)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+3)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+4)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+5)
    ixsecUL+=6
    myStructCmd += "}"
    rt.gROOT.ProcessLine(myStructCmd)
    from ROOT import MyStruct

    s = MyStruct()
    smoothXsecTree.Branch("mg", rt.AddressOf(s,"mg"),'mg/D')
    smoothXsecTree.Branch("mchi", rt.AddressOf(s,"mchi"),'mchi/D')

    ixsecUL = 0
    smoothXsecTree.Branch("xsecULObs_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+0)),'xsecUL%i/D'%(ixsecUL+0))
    smoothXsecTree.Branch("xsecULExpPlus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+1)),'xsecUL%i/D'%(ixsecUL+1))
    smoothXsecTree.Branch("xsecULExpPlus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+2)),'xsecUL%i/D'%(ixsecUL+2))
    smoothXsecTree.Branch("xsecULExp_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+3)),'xsecUL%i/D'%(ixsecUL+3))
    smoothXsecTree.Branch("xsecULExpMinus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+4)),'xsecUL%i/D'%(ixsecUL+4))
    smoothXsecTree.Branch("xsecULExpMinus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+5)),'xsecUL%i/D'%(ixsecUL+5))

    for mg in range(int(mgMin-12.5),int(mgMax+binWidth-12.5),int(binWidth)):
        for mchi in range(int(mchiMin-12.5),int(mchiMax+binWidth-12.5),int(binWidth)):
            s.mg = mg
            s.mchi = mchi
            
            xsecULObs = xsecUL["Obs"].GetBinContent(xsecUL["Obs"].FindBin(mg,mchi))
            xsecULExp = xsecUL["Exp"].GetBinContent(xsecUL["Exp"].FindBin(mg,mchi))
            xsecULExpPlus = xsecUL["ExpMinus"].GetBinContent(xsecUL["ExpMinus"].FindBin(mg,mchi))
            xsecULExpMinus = xsecUL["ExpPlus"].GetBinContent(xsecUL["ExpPlus"].FindBin(mg,mchi))
            xsecULExpPlus2 = xsecUL["ExpMinus2"].GetBinContent(xsecUL["ExpMinus2"].FindBin(mg,mchi))
            xsecULExpMinus2 = xsecUL["ExpPlus2"].GetBinContent(xsecUL["ExpPlus2"].FindBin(mg,mchi))
            
            exec 's.xsecUL%i = xsecULObs'%(ixsecUL+0)
            exec 's.xsecUL%i = xsecULExpPlus2'%(ixsecUL+1)
            exec 's.xsecUL%i = xsecULExpPlus'%(ixsecUL+2)
            exec 's.xsecUL%i = xsecULExp'%(ixsecUL+3)
            exec 's.xsecUL%i = xsecULExpMinus'%(ixsecUL+4)
            exec 's.xsecUL%i = xsecULExpMinus2'%(ixsecUL+5)

            if xsecULObs > 0.:
                smoothXsecTree.Fill()

    smoothOutFile.cd()
    smoothXsecTree.Write()
    smoothOutFile.Close()

    print len(gchipairs(model)), "total points"
    print len(toFix), "points to fix"
    print toFix
