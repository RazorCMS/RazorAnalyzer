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
from PlotFit import *


def binnedFit(pdf, data, fitRegion='Full'):
    #nll = pdf.createNLL(data)
    #m = rt.RooMinuit(nll)
    #m.migrad()
    #m.hesse()
    #fr = m.save()
    fr = extRazorPdf.fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Minimizer('Minuit2','improve'),rt.RooFit.Range(fitRegion),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.SumW2Error(False))
    #fr = extRazorPdf.fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Range(fitRegion),rt.RooFit.SumW2Error(False))
    
    return fr

def convertSideband(name,w,x,y,z):
    names = name.split(',')
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    iBinX = -1
    sidebandBins = []
    for ix in range(1,len(x)):
        for iy in range(1,len(y)):
            for iz in range(1,len(z)):
                iBinX+=1
                w.var('MR').setVal((x[ix]+x[ix-1])/2.)
                w.var('Rsq').setVal((y[iy]+y[iy-1])/2.)
                w.var('nBtag').setVal((z[iz]+z[iz-1])/2.)
                inSideband = 0
                for fitname in names:
                    inSideband += ( w.var('MR').inRange(fitname) * w.var('Rsq').inRange(fitname) * w.var('nBtag').inRange(fitname) )
                if inSideband: sidebandBins.append(iBinX)
                    
    sidebandGroups = []
    for k, g in groupby(enumerate(sidebandBins), lambda (i,x):i-x):
        consecutiveBins = map(itemgetter(1), g)
        sidebandGroups.append([consecutiveBins[0],consecutiveBins[-1]+1])

    newsidebands = ''
    minSideband = 0
    for iSideband in range(0,nBins):
        if w.var('th1x').hasRange('sideband%i'%iSideband):
            minSideband += 1
        
    for iSideband, sidebandGroup in enumerate(sidebandGroups):
        w.var('th1x').setRange("sideband%i"%(iSideband+minSideband),sidebandGroup[0],sidebandGroup[1])
        newsidebands+='sideband%i,'%(iSideband+minSideband)
    newsidebands = newsidebands[:-1]

    return newsidebands
 
    
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
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (useful for visualizing initial parameters)")
    parser.add_option('--fit-region',dest="fitRegion",default="Full",type="string",
                  help="Fit region")


    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    
    box = options.box
    lumi = options.lumi
    noFit = options.noFit
    fitRegion = options.fitRegion
    
    lumi_in = 0.

    data = None
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get('w'+box)
            data = workspace.data('RMRTree')
            lumi_in = 1000.*float([g.replace('lumi-','') for g in f.split('_') if g.find('lumi')!=-1][0])
    if data is None:
        print "give a root file as input"
  
    w = rt.RooWorkspace("w"+box)
    rootTools.Utils.importToWS(w,data)
    
    paramNames, bkgs = initializeWorkspace(w,cfg,box)
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    xFine = array('d', [x[0]+i*(x[-1]-x[0])/100. for i in range(0,101)]) # MR binning fine
    yFine = array('d', [y[0]+i*(y[-1]-y[0])/100. for i in range(0,101)]) # Rsq binning fine
    zFine = array('d', cfg.getBinning(box)[2]) # nBtag binning fine
    nBinsFine = (len(xFine)-1)*(len(yFine)-1)*(len(zFine)-1)

    
    th1x = w.var('th1x')
    sideband = convertSideband(fitRegion,w,x,y,z)
    plotRegion = 'Full'
    plotband = convertSideband(plotRegion,w,x,y,z)
    
    myTH1 = convertDataset2TH1(data, cfg, box, w)
    myTH1.Scale(lumi/lumi_in)
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), rt.RooFit.Import(myTH1))
    rootTools.Utils.importToWS(w,dataHist)
    
    setStyle()
    
    extRazorPdf = w.pdf('extRazorPdf')

    if noFit:
        fr = rt.RooFitResult()
    else:
        nll = extRazorPdf.createNLL(dataHist,rt.RooFit.Range(sideband))
        fr = binnedFit(extRazorPdf,dataHist,sideband)
        fr.Print('v')    
        rootTools.Utils.importToWS(w,fr)
        #rootTools.Utils.importToWS(w,nll)
        
    
    th1x.setBins(nBins)

    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())
    opt = []
    [opt.append(rt.RooFit.CutRange(myRange)) for myRange in plotband.split(',')]
    asimov_reduce = asimov.reduce(*opt)
    dataHist_reduce = dataHist.reduce(*opt)

    c = rt.TCanvas('c','c',500,400)
    
    rt.TH1D.SetDefaultSumw2()
    rt.TH2D.SetDefaultSumw2()
    rt.TH3D.SetDefaultSumw2()

    h_th1x = asimov_reduce.createHistogram('h_th1x',th1x)
    h_data_th1x = dataHist_reduce.createHistogram('h_data_th1x',th1x)
    
    h_data_nBtagRsqMR = get3DHistoFrom1D(h_data_th1x,x,y,z,"h_data_nBtagRsqMR")
    h_nBtagRsqMR = get3DHistoFrom1D(h_th1x,x,y,z,"h_nBtagRsqMR")
                
    h_data_RsqMR = h_data_nBtagRsqMR.Project3D("yxe")
    h_data_MR = h_data_nBtagRsqMR.Project3D("xe")
    h_data_Rsq = h_data_nBtagRsqMR.Project3D("ye")
    h_RsqMR = h_nBtagRsqMR.Project3D("yxe")
    h_MR = h_nBtagRsqMR.Project3D("xe")
    h_Rsq = h_nBtagRsqMR.Project3D("ye")

    
    if len(z)>1:
        h_MR_components = []
        h_Rsq_components = []
        h_labels = []        
        h_colors = [rt.kOrange,rt.kViolet,rt.kRed,rt.kGreen]
        for k in range(1,len(z)):
            h_MR_components.append(h_nBtagRsqMR.ProjectionX("MR_%ibtag"%z[k-1],0,-1,k,k,""))
            h_Rsq_components.append(h_nBtagRsqMR.ProjectionY("Rsq_%ibtag"%z[k-1],0,-1,k,k,""))
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

    
         
    print1DProj(c,h_MR,h_data_MR,options.outDir+"/h_MR_%s.pdf"%box,"M_{R} [GeV]","Events",lumiLabel,boxLabel,None,h_MR_components,h_colors,h_labels)
    print1DProj(c,h_Rsq,h_data_Rsq,options.outDir+"/h_Rsq_%s.pdf"%box,"R^{2}","Events",lumiLabel,boxLabel,None,h_Rsq_components,h_colors,h_labels)

    h_data_nBtagRsqMR_fine = rt.TH3D("h_data_nBtagRsqMR_fine","h_data_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    h_nBtagRsqMR_fine = rt.TH3D("h_nBtagRsqMR_fine","h_nBtagRsqMR_fine",len(xFine)-1,xFine,len(yFine)-1,yFine,len(zFine)-1,zFine)
    w.data("RMRTree").fillHistogram(h_data_nBtagRsqMR_fine,rt.RooArgList(w.var("MR"),w.var("Rsq"),w.var("nBtag")))
    for i in range(1,len(xFine)):
        for j in range(1,len(yFine)):
            for k in range(1,len(zFine)):
                value, errorFlag= getBinEvents(i,j,k,xFine,yFine,zFine,w,box)
                if not errorFlag:
                    h_nBtagRsqMR_fine.SetBinContent(i,j,k,value)

            
    h_data_RsqMR_fine = h_data_nBtagRsqMR_fine.Project3D("yxe")
    h_RsqMR_fine = h_nBtagRsqMR_fine.Project3D("yxe")
    
    
    print2DScatter(c,h_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Fit",lumiLabel,boxLabel,x,y,h_data_RsqMR_fine.GetMinimum(),h_data_RsqMR_fine.GetMaximum())
    print2DScatter(c,h_data_RsqMR_fine,options.outDir+"/h_RsqMR_scatter_data_log_%s.pdf"%(box),"M_{R} [GeV]", "R^{2}", "Sim. Data",lumiLabel,boxLabel,x,y,h_data_RsqMR_fine.GetMinimum(),h_data_RsqMR_fine.GetMaximum())
    

    outFileName = "BinnedFitResults_%s.root"%(box)
    outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
