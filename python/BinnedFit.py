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
import RunToys


def binnedFit(pdf, data):
    #nll = pdf.createNLL(data)
    #m = rt.RooMinuit(nll)
    #m.migrad()
    #m.hesse()
    #fr = m.save()
    fr = extRazorPdf.fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Minimizer('Minuit2','improve'),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.SumW2Error(False))
    
    return fr
 
def setStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store cards")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    
    box = options.box
    lumi = options.lumi
    
    lumi_in = 0.
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get('w'+box)
            data = workspace.data('RMRTree')
            lumi_in = 1000.*float([g.replace('lumi-','') for g in f.split('_') if g.find('lumi')!=-1][0])

  
    w = rt.RooWorkspace("w"+box)
    
    paramNames = initializeWorkspace(w,cfg,box)
    rootTools.Utils.importToWS(w,data)
    
    th1x = w.var('th1x')
    
    myTH1 = convertDataset2TH1(data, cfg, box, w)
    myTH1.Scale(lumi/lumi_in)
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), myTH1)
    rootTools.Utils.importToWS(w,dataHist)

    setStyle()
    
    extRazorPdf = w.pdf('extRazorPdf')

    fr = binnedFit(extRazorPdf,dataHist)
    #fr = rt.RooFitResult()
    fr.Print('v')
    
    rootTools.Utils.importToWS(w,fr)
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    th1x.setBins(nBins)

    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())

    c = rt.TCanvas('c','c',500,400)
    
    rt.TH1D.SetDefaultSumw2()
    rt.TH2D.SetDefaultSumw2()
    rt.TH3D.SetDefaultSumw2()

    h_th1x = asimov.createHistogram('h_th1x',th1x)    
    h_nBtagRsqMR = rt.TH3D("h_nBtagRsqMR","h_nBtagRsqMR",len(x)-1,x,len(y)-1,y,len(z)-1,z)

    h_data_th1x = dataHist.createHistogram('h_data_th1x',th1x)
    h_data_nBtagRsqMR = rt.TH3D("h_data_nBtagRsqMR","h_data_nBtagRsqMR",len(x)-1,x,len(y)-1,y,len(z)-1,z)

    iBinX = -1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                h_data_nBtagRsqMR.SetBinContent(i,j,k,h_data_th1x.GetBinContent(iBinX+1))
                h_data_nBtagRsqMR.SetBinError(i,j,k,h_data_th1x.GetBinError(iBinX+1))
                h_nBtagRsqMR.SetBinContent(i,j,k,h_th1x.GetBinContent(iBinX+1))
                h_nBtagRsqMR.SetBinError(i,j,k,h_th1x.GetBinError(iBinX+1))

                
    h_data_RsqMR = h_data_nBtagRsqMR.Project3D("yxe")
    h_data_MR = h_data_nBtagRsqMR.Project3D("xe")
    h_data_Rsq = h_data_nBtagRsqMR.Project3D("ye")
    h_RsqMR = h_nBtagRsqMR.Project3D("yxe")
    h_MR = h_nBtagRsqMR.Project3D("xe")
    h_Rsq = h_nBtagRsqMR.Project3D("ye")

    


    btagLabel = ""
    if z[-1] == z[0]+1 and z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    elif z[-1] == z[0]+1:
        btagLabel = "%i b-tag" % z[0]
    elif z[-1]==4:
        btagLabel = "#geq %i b-tag" % z[0]
    else:
        btagLabel = "%i-%i b-tag" % (z[0],z[-2])

    lumiLabel = "%i fb^{-1} (13 TeV)" % (lumi/1000.)
    boxLabel = "razor %s %s" % (box,btagLabel)
         
    RunToys.print1DCanvas(c,h_MR,h_data_MR,options.outDir+"/h_MR.pdf","M_{R} [GeV]","Events",lumiLabel,boxLabel)
    RunToys.print1DCanvas(c,h_Rsq,h_data_Rsq,options.outDir+"/h_Rsq.pdf","R^{2}","Events",lumiLabel,boxLabel)
    

    outFileName = "BinnedFitResults_%s.root"%(box)
    outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
