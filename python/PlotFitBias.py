from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
from itertools import *
from operator import *
from WriteDataCard import initializeWorkspace, convertDataset2TH1
import os
import random
import sys
import math
from PlotFit import setStyle,print1DProj,print2DScatter,get3DHistoFrom1D,getBinEvents,convertSideband,densityCorr,getBinSumDicts

def getTrueValue(varName,toyTree):
    
    toyTree.GetEntry(0)
    trueValue = eval(varName.replace('b','toyTree.b'))
    
    return trueValue

def getBiasHistos(varName,toyTree,name):
    
    toyTree.Draw('%s>>htest%s'%(varName,name),'migrad_ff==0&&migrad_sf==0&&toy_num>-1')
    #toyTree.Draw('%s>>htest%s'%(varName,name),'toy_num>-1')
    htemp = rt.gPad.GetPrimitive("htest%s"%name)
    #xmax = htemp.GetMean()+3.*htemp.GetRMS()
    #xmin = htemp.GetMean()-3.*htemp.GetRMS()
    if htemp.GetXaxis().GetXmin()>=0:
        xmin = htemp.GetXaxis().GetXmin()
        xmax = htemp.GetXaxis().GetXmax()
    else:        
        #xmax = max(abs(htemp.GetXaxis().GetXmax()),abs(htemp.GetXaxis().GetXmin()))
        #xmin = -xmax        
        xmin = max(htemp.GetXaxis().GetXmin(),htemp.GetMean()-3.*htemp.GetRMS())
        xmax = min(htemp.GetXaxis().GetXmax(),htemp.GetMean()+3.*htemp.GetRMS())
        
    
    
    h = rt.TH1D('h_%s'%name,'h_%s'%name,20,xmin,xmax)
    toyTree.Project('h_%s'%name,varName)

    return h

def printBiasHisto(h_data,xTitle,yTitle,lumi,boxLabel,btagLabel,printName):
        h_data.SetMarkerColor(rt.kBlack)
        h_data.SetMarkerStyle(20)
        h_data.SetLineColor(rt.kBlack)
        h_data.GetXaxis().SetTitle(xTitle)
        h_data.GetYaxis().SetTitle(yTitle)
        h_data.GetXaxis().SetLabelSize(0.05)
        h_data.GetYaxis().SetLabelSize(0.05)
        h_data.GetXaxis().SetTitleSize(0.06)
        h_data.GetYaxis().SetTitleSize(0.06)
        h_data.GetXaxis().SetTitleOffset(1)
        h_data.GetYaxis().SetTitleOffset(1)
        h_data.SetMaximum(1.7*h_data.GetBinContent(h_data.GetMaximumBin()))
            
        h_data.Draw("pe")

        tLeg = rt.TLegend(0.65,0.55,0.89,0.75)
        tLeg.SetFillColor(0)
        tLeg.SetTextFont(42)
        tLeg.SetLineColor(0)
        tLeg.AddEntry(h_data,"Toy Datasets","lep")
        tLeg.AddEntry(None,"Mean = %.2f"%h_data.GetMean(),"")
        tLeg.AddEntry(None,"RMS  = %.2f"%h_data.GetRMS(),"")
            
        tLeg.Draw("same")

        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.05)
        l.SetTextFont(42)
        l.SetNDC()
        l.DrawLatex(0.15,0.9,"CMS private")
        
        lumiLabel = "%.1f fb^{-1} (13 TeV)" % (lumi/1000)
        l.DrawLatex(0.7,0.9,"%s"%lumiLabel)
        l.SetTextFont(52)
        l.SetTextSize(0.045)
        l.DrawLatex(0.2,0.82,boxLabel)
        l.DrawLatex(0.3,0.77,btagLabel)
        c.Print(printName)

            
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
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
                  help="input fit file")    
    parser.add_option('-t','--input-toy-file',dest="inputToyFile", default=None,type="string",
                  help="input toy file")
    
    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)    
    lumi = options.lumi
    inputFitFile = options.inputFitFile
    box = options.box

    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)

    

    ixMin = 3
    iyMin = 3
    if box in ['MuMultiJet','EleMultiJet']:
        ixMin = 2
        iyMin = 3
    binSumDict = getBinSumDicts("z", ixMin, len(x)-1, iyMin, len(y)-1, 0, len(z)-1, x, y, z)

    
    toyTree = None
    if options.inputToyFile is not None:
        toyFiles = options.inputToyFile.split(',')
        toyTree = rt.TChain("myTree")
        for toyFile in toyFiles:
            toyTree.Add(toyFile)

    biasHistos = []

    setStyle()
    c = rt.TCanvas('c','c',500,400)
    
    c.SetLeftMargin(0.15)
    c.SetRightMargin(0.1)
    c.SetTopMargin(0.12)
    c.SetBottomMargin(0.15)
    
    boxLabel = 'razor %s box'%box
    for iz in range(1,len(z)):
        btagLabel = 'M_{R} > %i GeV, R^{2} > %.2f, '%(x[ixMin-1],y[iyMin-1])
        if z[iz]==4:
            btagLabel += "#geq %i b-tag" % z[iz-1]
        else:            
            btagLabel += "%i b-tag" % z[iz-1]
            
        sfBinSum =  binSumDict[iz].replace('+','_sf+')+'_sf'
        ffBinSum =  binSumDict[iz].replace('+','_ff+')+'_ff'
        toyBinSum =  binSumDict[iz].replace('+','_toy+')+'_toy'
        trueValue = getTrueValue(toyBinSum,toyTree)
        h_sf_ff = getBiasHistos(sfBinSum+'-('+ffBinSum+')',toyTree,'h_sf_ff_%ibtag'%z[iz-1])
        h_sf = getBiasHistos(sfBinSum,toyTree,'h_sf_%ibtag'%z[iz-1])
        h_ff = getBiasHistos(ffBinSum,toyTree,'h_ff_%ibtag'%z[iz-1])
        h_toy = getBiasHistos(toyBinSum,toyTree,'h_toy_%ibtag'%z[iz-1])
        h_sf_ff_divtoy = getBiasHistos('('+sfBinSum+'-('+ffBinSum+'))/('+toyBinSum+')',toyTree,'h_sf_ff_divtoy_%ibtag'%z[iz-1])
        h_sf_ff_divff = getBiasHistos('('+sfBinSum+'-('+ffBinSum+'))/('+ffBinSum+')',toyTree,'h_sf_ff_divff_%ibtag'%z[iz-1])
        h_sf_ff_divtrue = getBiasHistos('('+sfBinSum+'-('+ffBinSum+'))/(%f)'%trueValue,toyTree,'h_sf_ff_divtrue_%ibtag'%z[iz-1])

        printBiasHisto(h_sf,'Sideband Fit Yield','Toy Datasets',lumi,boxLabel,btagLabel,options.outDir+"/yield_sf_%ibtag_%s.pdf"%(z[iz-1],box))
        printBiasHisto(h_ff,'Full Fit Yield','Toy Datasets',lumi,boxLabel,btagLabel,options.outDir+"/yield_ff_%ibtag_%s.pdf"%(z[iz-1],box))
        printBiasHisto(h_toy,'Toy Yield','Toy Datasets',lumi,boxLabel,btagLabel,options.outDir+"/yield_toy_%ibtag_%s.pdf"%(z[iz-1],box))
        printBiasHisto(h_sf_ff,'Sideband Fit - Full Fit','Toy Datasets',lumi,boxLabel,btagLabel,options.outDir+"/bias_sf_ff_%ibtag_%s.pdf"%(z[iz-1],box))
        printBiasHisto(h_sf_ff_divtoy,'(Sideband Fit - Full Fit) / Toy','Toy Datasets',lumi,boxLabel,btagLabel,options.outDir+"/bias_sf_ff_divtoy_%ibtag_%s.pdf"%(z[iz-1],box))
        printBiasHisto(h_sf_ff_divff,'(Sideband Fit - Full Fit) / Full Fit','Toy Datasets',lumi,boxLabel,btagLabel,options.outDir+"/bias_sf_ff_divff_%ibtag_%s.pdf"%(z[iz-1],box))
        printBiasHisto(h_sf_ff_divtrue,'(Sideband Fit - Full Fit) / True Value','Toy Datasets',lumi,boxLabel,btagLabel,options.outDir+"/bias_sf_ff_divtrue_%ibtag_%s.pdf"%(z[iz-1],box))
