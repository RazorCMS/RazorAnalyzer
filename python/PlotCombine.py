from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="box", default="T1bbbb",type="string",
                  help="signal model name")
    parser.add_option('--mGluino',dest="mGluino", default=-1,type="float",
                  help="mass of gluino")
    parser.add_option('--mStop',dest="mStop", default=-1,type="float",
                  help="mass of stop")
    parser.add_option('--mLSP',dest="mLSP", default=-1,type="float",
                  help="mass of LSP")
    parser.add_option('--lumi-array',dest="lumi_array", default="4",type="string",
                  help="lumi array in fb^-1, e.g.: 0.2,4,10")

    (options,args) = parser.parse_args()


    lumiArray = array('d',[float(lumi) for lumi in options.lumi_array.split(',')])
    box = options.box
    model = options.model
    mGluino = options.mGluino
    mLSP = options.mLSP
    mStop = options.mStop
    
    if mGluino!=-1:
        massPoint = "mGl-%i_mLSP-%i"%(mGluino,mLSP)
    if mStop!=-1:
        massPoint = "mStop-%i_mLSP-%i"%(mStop,mLSP)
    
    expArray = array('d')
    expp1sigmaArray = array('d')
    expm1sigmaArray = array('d')
    expp2sigmaArray = array('d')
    expm2sigmaArray = array('d')
    zeroArray = array('d',[0 for lumi in lumiArray])
    xsecArray = array('d',[thyXsec for lumi in lumiArray])
    xsecp1sigmaArray = array('d',[thyXsec*(thyXsecErr)/2. for lumi in lumiArray])
    xsecm1sigmaArray = array('d',[thyXsec*(thyXsecErr)/2. for lumi in lumiArray])
    
    for lumi in lumiArray:
        tfile = rt.TFile.Open('cards/higgsCombine%s_%s_lumi-%s_%s.Asymptotic.mH120.root'%(model,massPoint,lumi,box))
        limit = tfile.Get('limit')
        limit.Draw('>>elist','','entrylist')
        elist = rt.gDirectory.Get('elist')
        entry = elist.Next()
        limits = []
        while entry>=0:
            limit.GetEntry(entry)
            limits.append(limit.limit)
            entry = elist.Next()
        expArray.append(limits[2]*refXsec)
        expp2sigmaArray.append((limits[4]-limits[2])*refXsec)
        expp1sigmaArray.append((limits[3]-limits[2])*refXsec)
        expm1sigmaArray.append((limits[2]-limits[1])*refXsec)
        expm2sigmaArray.append((limits[2]-limits[0])*refXsec)

    print expp1sigmaArray
    print expm1sigmaArray
    print expp2sigmaArray
    print expm2sigmaArray
    
    
    expGraph = rt.TGraph(len(lumiArray),lumiArray,expArray)
    expGraph.SetLineStyle(2)
    
    exp2sigmaGraph = rt.TGraphAsymmErrors(len(lumiArray), lumiArray, expArray ,zeroArray,zeroArray, expm2sigmaArray,expp2sigmaArray)
    exp2sigmaGraph.SetLineColor(5)
    exp2sigmaGraph.SetFillColor(5)
    exp2sigmaGraph.SetFillStyle(1001)

    exp1sigmaGraph = rt.TGraphAsymmErrors(len(lumiArray), lumiArray, expArray ,zeroArray,zeroArray, expm1sigmaArray,expp1sigmaArray)
    exp1sigmaGraph.SetLineColor(rt.kGreen-7)
    exp1sigmaGraph.SetFillColor(rt.kGreen-7)
    exp1sigmaGraph.SetFillStyle(1001)

    xsecGraphErr = rt.TGraphAsymmErrors(len(lumiArray),lumiArray,xsecArray,zeroArray,zeroArray,xsecm1sigmaArray,xsecp1sigmaArray)
    xsecGraphErr.SetFillStyle(1001)
    xsecGraphErr.SetLineColor(rt.kOrange)
    xsecGraphErr.SetFillColor(rt.kBlue-7)
    
    xsecGraph = rt.TGraph(len(lumiArray),lumiArray,xsecArray)
    xsecGraph.SetMarkerSize(0)
    xsecGraph.SetLineStyle(1)
    xsecGraph.SetLineColor(rt.kOrange)

    
    h_limit = rt.TMultiGraph()
    c = rt.TCanvas('c','c',600,400)
    c.SetLogy()
    h_limit.Add(exp2sigmaGraph)
    h_limit.Add(exp1sigmaGraph)
    h_limit.Add(xsecGraphErr)
    h_limit.Draw("a3")
    h_limit.GetXaxis().SetLimits(0.2,5.)
    h_limit.GetXaxis().SetTitle("Integrated Luminosity [fb^{-1}]")
    h_limit.GetYaxis().SetTitle("95% C.L. Upper Limit Cross Section [fb]")
    h_limit.SetMaximum(1.e3)
    #h_limit.SetMaximum(1.e4)
    h_limit.SetMinimum(2.)
    #h_limit.SetMinimum(50.)
    h_limit.Draw("a3")
    
    expGraph.Draw("c same")
    xsecGraph.Draw("c same")
    

    demo = exp1sigmaGraph.Clone()
    demo.SetLineColor(rt.kBlack)
    demo.SetLineStyle(2)

    
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.045)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.18,0.84,"CMS Simulation (13 TeV)")
    l.DrawLatex(0.18,0.77,"Razor %s Box"%box)
    if model=="T1bbbb":
        l.DrawLatex(0.53,0.84,"pp #rightarrow #tilde{g}#tilde{g},  #tilde{g}#rightarrowb#bar{b}#tilde{#chi}^{0}_{1}")
        l.DrawLatex(0.53,0.77,"m_{#tilde{g}} = %i GeV, m  _{#tilde{#chi}} = %i GeV"%(mGluino,mLSP))
        
    leg = rt.TLegend(0.52,0.54,0.89,0.74)
    leg.SetTextFont(42)
    leg.SetFillColor(rt.kWhite)
    leg.SetLineColor(rt.kWhite)
    leg.AddEntry(xsecGraphErr, "#sigma_{NLO+NLL} (#tilde{g}#tilde{g}) #pm 1 #sigma (theory)","lf")
    leg.AddEntry(demo, "expected #pm 1 #sigma (experiment)","lf")

    leg.Draw()
    
    if model.find("T1")!=-1:
        c.Print("%s/limit_%s_%i_%i.pdf"%(outDir,model,mGluino,mLSP))
        c.Print("%s/limit_%s_%i_%i.C"%(outDir,model,mGluino,mLSP))
    elif model.find("T2")!=-1:
        c.Print("%s/limit_%s_%i_%i.pdf"%(outDir,model,mStop,mLSP))
        c.Print("%s/limit_%s_%i_%i.C"%(outDir,model,mStop,mLSP))
