import ROOT as rt
import sys
import rootTools
import glob
from math import *
import os
from array import *
import random
from optparse import OptionParser
from framework import Config
from PlotFit import *


def testWorkspace(w,outFile,singleBox,options):
     
    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
    
    w.Print("v")
    CMS_th1x = w.var("th1x")
    CMS_channel = w.cat("CMS_channel")
    
    CMS_set = rt.RooArgSet()
    CMS_set.add(CMS_channel)
    CMS_set.add(CMS_th1x)
    
    data_obs = w.data("data_obs")

    #model_b = w.pdf("model_b")
    model_b = w.pdf("pdf_bin%s_bonly"%singleBox)
    #model_s = w.pdf("model_s")
    model_s = w.pdf("pdf_bin%s"%singleBox)

    #data_obs = model_b.generateBinned(CMS_set,rt.RooFit.Asimov())
    
    signalNuis = ["lumi"]
            
    for varName in signalNuis:
        print varName
        #w.var(varName).setConstant()

    r = w.var('r')
    r.setMax(20.)
    poi = w.set('POI')
    
    allParams = model_b.getParameters(data_obs)
    rt.RooStats.RemoveConstantParameters(allParams)

    opt = rt.RooLinkedList()
    opt.Add(rt.RooFit.CloneData(False))
    opt.Add(rt.RooFit.Constrain(allParams))

    nll = model_s.createNLL(data_obs,opt)
    
    minim = rt.RooMinimizer(nll)

    strategy = rt.Math.MinimizerOptions.DefaultStrategy()

    minim.setStrategy(strategy)
    tol = rt.Math.MinimizerOptions.DefaultTolerance()
    #tol = max(tol,1.0)
    tol = 0.001
    minim.setEps(tol)

    minim.optimizeConst(2)

    minimizer = rt.Math.MinimizerOptions.DefaultMinimizerType()
    algorithm = rt.Math.MinimizerOptions.DefaultMinimizerAlgo()

    status = minim.minimize(minimizer,algorithm)
    
    fr_s = model_s.fitTo(data_obs,rt.RooFit.Save())
    fr_s.Print("v")
    nll_s = fr_s.minNll()

    rBestFit = r.getVal()
    r.setVal(0.)

    asimov_b = model_s.generateBinned(CMS_set,rt.RooFit.Asimov())
    h_bkgd_th1x = asimov_b.createHistogram('h_th1x',CMS_th1x)
    
    r.setVal(rBestFit)
    #r.setVal(0.742214)
    #r.setConstant()    
    #fr_s = model_s.fitTo(data_obs,rt.RooFit.Save())
    
    norm = {}
    for bkg in [0, 1, 2, 3]:
        if w.var('shapeBkg_%s_TTj%ib_%s__norm'%(singleBox,bkg,singleBox))!=None:            
            norm[bkg] = w.var('shapeBkg_%s_TTj%ib_%s__norm'%(singleBox,bkg,singleBox)).getVal()
            w.var('shapeBkg_%s_TTj%ib_%s__norm'%(singleBox,bkg,singleBox)).setVal(0.)

            
    asimov_s = model_s.generateBinned(CMS_set,rt.RooFit.Asimov())
    h_sig_th1x = asimov_s.createHistogram('h_th1x',CMS_th1x)
    
    for bkg in [0, 1, 2, 3]:
        if w.var('shapeBkg_%s_TTj%ib_%s__norm'%(singleBox,bkg,singleBox))!=None:            
            w.var('shapeBkg_%s_TTj%ib_%s__norm'%(singleBox,bkg,singleBox)).setVal(norm[bkg])            
            
    cfg = Config.Config(options.config)    
    x = array('d', cfg.getBinning(singleBox)[0]) # MR binning
    y = array('d', cfg.getBinning(singleBox)[1]) # Rsq binning
    z = array('d', cfg.getBinning(singleBox)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    h_data_th1x = data_obs.createHistogram('h_data_th1x',CMS_th1x)
    
    h_data_nBtagRsqMR = get3DHistoFrom1D(h_data_th1x,x,y,z,"h_data_nBtagRsqMR")
    h_bkgd_nBtagRsqMR = get3DHistoFrom1D(h_bkgd_th1x,x,y,z,"h_bkgd_nBtagRsqMR")
    h_sig_nBtagRsqMR = get3DHistoFrom1D(h_sig_th1x,x,y,z,"h_sig_nBtagRsqMR")
    
    h_data_MR = h_data_nBtagRsqMR.Project3D("xe")
    h_data_Rsq = h_data_nBtagRsqMR.Project3D("ye")
    h_bkgd_MR = h_bkgd_nBtagRsqMR.Project3D("xe")
    h_bkgd_Rsq = h_bkgd_nBtagRsqMR.Project3D("ye")
    h_sig_MR = h_sig_nBtagRsqMR.Project3D("xe")
    h_sig_Rsq = h_sig_nBtagRsqMR.Project3D("ye")

    h_total_MR = h_bkgd_MR.Clone("h_total_MR")
    h_total_MR.Add(h_sig_MR)
    h_colors = [rt.kGreen,rt.kRed]
    h_MR_components = [h_bkgd_MR, h_sig_MR]
    h_labels = ['bkgd','signal']

    
    h_total_Rsq = h_bkgd_Rsq.Clone("h_total_Rsq")
    h_total_Rsq.Add(h_sig_Rsq)
    h_colors = [rt.kGreen,rt.kRed]
    h_Rsq_components = [h_bkgd_Rsq, h_sig_Rsq]
    h_labels = ['bkgd','signal']


    btagLabel = "#geq %i b-tag" % z[0]
    lumiLabel = "%.0f pb^{-1} (13 TeV)" % (options.lumi)
    boxLabel = "razor %s %s" % (singleBox,btagLabel)
    dataString = "Data"
    
    c = rt.TCanvas('d','d',500,400)  
    print1DProj(c,outFile,h_total_MR,h_data_MR,"./h_MR_%s.pdf"%singleBox,"M_{R} [GeV]","Events",lumiLabel,boxLabel,True,None,h_MR_components,h_colors,h_labels)
    print1DProj(c,outFile,h_total_Rsq,h_data_Rsq,"./h_Rsq_%s.pdf"%singleBox,"R^{2}","Events",lumiLabel,boxLabel,True,None,h_Rsq_components,h_colors,h_labels)
    
    pll = nll.createProfile(poi)
    n2ll = rt.RooFormulaVar("n2ll","2*@0-2*%f"%nll_s,rt.RooArgList(nll))
    p2ll = n2ll.createProfile(poi)
    
    print "signal+background nll = %f at r = %f"%(nll_s,r.getVal())
    #minim.optimizeConst(False)

    
    #c = rt.TCanvas('c','c',500,500)
    frame = r.frame(rt.RooFit.Bins(10),rt.RooFit.Range(0.0,options.rMax),rt.RooFit.Title("r frame"))
    frame.SetMinimum(0)
    frame.SetMaximum(6)

    
    #plot_opt = rt.RooLinkedList()
    #plot_opt.Add(rt.RooFit.ShiftToZero())
    #plot_opt.Add(rt.RooFit.LineColor(rt.kBlack))
    
    n2ll.plotOn(frame,rt.RooFit.ShiftToZero(),rt.RooFit.LineStyle(2),rt.RooFit.Name("n2ll"))
    #p2ll.plotOn(frame,rt.RooFit.LineColor(rt.kBlack),rt.RooFit.Name("p2ll"))

    tlines = []
    cl = 0.95
    crossing = rt.TMath.Power(rt.Math.normal_quantile(1-0.5*(1-cl), 1.0),2)
    tline = rt.TLine(0,crossing,options.rMax,crossing)
    tline.SetLineColor(rt.kRed)
    tline.SetLineWidth(2)
    tlines.append(tline)
    for tline in tlines:
        frame.addObject(tline,"")

        
    frame.Draw()
    frame.SetMinimum(0)
    frame.SetMaximum(6)
    
    frame.SetXTitle("r (signal strength)")
    frame.SetYTitle("-2 #Delta log L(asimov)")
    frame.SetTitleSize(0.04,"X")
    frame.SetTitleOffset(0.85,"X")
    frame.SetTitleSize(0.04,"Y")
    frame.SetTitleOffset(0.8,"Y")
    frame.SetLabelSize(0.04,"X")
    frame.SetLabelSize(0.04,"Y")
    frame.SetNdivisions(505,"X")
    
    leg = rt.TLegend(0.7,0.15,0.89,0.3)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)
    leg.AddEntry("p2ll", "stat + syst","l")
    leg.AddEntry("n2ll", "stat only","l")
    leg.Draw()
    
    l = rt.TLatex()
    l.SetTextAlign(12)
    l.SetTextSize(0.04)
    l.SetTextFont(42)
    l.SetNDC()
    l.SetTextSize(0.04)
    l.DrawLatex(0.12,0.85,"CMS Simulation, %.1f fb^{-1} (13 TeV)"%(options.lumi/1000))
    l.DrawLatex(0.12,0.80, "razor %s"%singleBox)    
    if options.mGluino!=-1:
        for line in open('data/gluino13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(options.mGluino))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1])*1000 # convert pb to fb
    if options.mStop!=-1:
        for line in open('data/stop13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(options.mStop))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1])*1000 # convert pb to fb
    if options.mGluino!=-1:
        if options.model=="T1tttt":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow t#bar{t}#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
        if options.model=="T1bbbb":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow b#bar{b}#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
        if options.model=="T1qqqq":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{g}#tilde{g}, #tilde{g}#rightarrow q#bar{q}#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
        l.DrawLatex(0.12,0.69,"m_{#tilde{g}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mGluino,options.mLSP))
    if options.mStop!=-1:
        if options.model=="T2tt":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{t}#tilde{t}, #tilde{t}#rightarrow t#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
            l.DrawLatex(0.12,0.69,"m_{#tilde{t}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mStop,options.mLSP))
        if options.model=="T2bb":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{b}#tilde{b}, #tilde{b}#rightarrow b#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
            l.DrawLatex(0.12,0.69,"m_{#tilde{b}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mStop,options.mLSP))
        if options.model=="T2qq":
            l.DrawLatex(0.12,0.75, "pp#rightarrow#tilde{q}#tilde{q}, #tilde{q}#rightarrow q#tilde{#chi}^{0}, #sigma = %.1f fb"%thyXsec)
            l.DrawLatex(0.12,0.69,"m_{#tilde{q}} = %i GeV, m_{#tilde{#chi}} = %i GeV"%(options.mStop,options.mLSP))
    outFile.cd()
    c.Write()
    frame.Write()
    outFile.Close()
    

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1tttt",type="string",
                  help="signal model name")
    parser.add_option('--mGluino',dest="mGluino", default=-1,type="float",
                  help="mass of gluino")
    parser.add_option('--mStop',dest="mStop", default=-1,type="float",
                  help="mass of stop")
    parser.add_option('--mLSP',dest="mLSP", default=-1,type="float",
                  help="mass of LSP")
    parser.add_option('--rMax',dest="rMax", default=3,type="float",
                  help="maximum of r (signal strength) in profile likelihood plot")
    parser.add_option('-o','--output',dest="output",default="./test.root",type="string",
                  help="Output file to store results")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")

    (options,args) = parser.parse_args()

     
    rt.gSystem.Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so")
    
    for f in args:
        if f.lower().endswith('.root'):
            workspaceFileName = f
    workspaceFile = rt.TFile.Open(workspaceFileName,"READ")
    outFile = rt.TFile.Open(options.output,"RECREATE")
    w = workspaceFile.Get("w")
    singleBox = options.box

    print "SINGLE BOX %s"%singleBox 
    testWorkspace(w,outFile,singleBox,options)
