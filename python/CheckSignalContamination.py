from optparse import OptionParser
import ROOT as rt
import sys
from array import *

#local imports
import rootTools
from framework import Config
from DustinTuple2RooDataSet import initializeWorkspace
from DustinTuples2DataCard import convertTree2TH1, uncorrelate
from RunCombine import exec_me
from macro.razorAnalysis import xbinsSignal, colsSignal
from macro.macro import importHists, makeTH2PolyFromColumns, fillTH2PolyFromTH2, stitch
import os

SIGNAL_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForMoriond20160124/combined"
#SIGNAL_DIR = "Signals"
BACKGROUND_DIR = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorMADD2015"
#BACKGROUND_DIR = "Moriond2016"

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--lumi-in',dest="lumi_in", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('--merge-bins',dest="mergeBins", action="store_true",
                  help="merge some bins in Rsq")    
    parser.add_option('-m','--model', default="T1bbbb", help="signal model name")
    parser.add_option('--mGluino',default=-1,type=int, help="mass of gluino")
    parser.add_option('--mLSP',default=-1,type=int, help="mass of LSP")
    (options,args) = parser.parse_args()


    cfg = Config.Config(options.config)

    box = options.box
    treeName = 'RazorInclusive'
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    unrollBins = None
    if options.mergeBins:
        btagBins = cfg.getBinning(box)[2][:-1]
        unrollBins = [(xbinsSignal[box][str(int(btags))+'B'], colsSignal[box][str(int(btags))+'B']) for btags in btagBins]

    mergeBinsString = ''
    if options.mergeBins:
        mergeBinsString = '--merge-bins'
        
    os.system('python python/SMSTemplates.py %s -c %s -d %s/ --lumi %s --box %s --no-signal-sys %s/SMS-%s_%i_%i.root'%(mergeBinsString,options.config,options.outDir,options.lumi,options.box,SIGNAL_DIR,options.model,options.mGluino,options.mLSP))
        
    signalFile = rt.TFile.Open('%s/SMS-%s_%i_%i_lumi-%.3f_%i-%ibtag_%s.root'%(options.outDir,options.model,options.mGluino,options.mLSP,options.lumi*1./1000.,z[0],z[-1]-1,box))
    sigTH1 = signalFile.Get('%s_%s'%(box,options.model))
    
    bkgdHistDict = importHists('%s/controlHistograms%s.root'%(BACKGROUND_DIR,box.replace('ControlRegion','').replace('WJet','WJetsSingleLepton')))
    tempTH2 = bkgdHistDict['Data'][('MR','Rsq')]

    
    #unroll into TH1F
    if unrollBins is None:
        nBins = (len(x)-1)*(len(y)-1)
        maxBins = nBins
        myTH1 = rt.TH1F(treeName,treeName,maxBins,0,maxBins)
        myTH1.SetDirectory(0) #prevent it from going out of scope
        myTH1.Sumw2()
        i = 0
        for ix in range(1,len(x)):
            for iy in range(1,len(y)):
                i += 1
                myTH1.SetBinContent(i,tempTH2.GetBinContent(ix,iy))
                myTH1.SetBinError(i,tempTH2.GetBinError(ix,iy))
    else:        
        print "Merging bins according to custom (MR-dependent) binning"
        layers = []
        #turn it into a TH2Poly with the reduced binning
        unrollRows = unrollBins[0][0]
        unrollCols = unrollBins[0][1]
        poly = makeTH2PolyFromColumns(tempTH2.GetName()+"poly", 'poly', unrollRows, unrollCols)
        fillTH2PolyFromTH2(tempTH2, poly)
        numbins = poly.GetNumberOfBins()
        unrolledSlice = rt.TH1D(tempTH2.GetName()+"Unroll", "slice", numbins, 0, numbins)
        for bn in range(1, numbins+1):
            unrolledSlice.SetBinContent(bn, poly.GetBinContent(bn))
            unrolledSlice.SetBinError(bn, poly.GetBinError(bn))
        layers.append(unrolledSlice)
        poly.Delete()
        myTH1 = stitch(layers)
        myTH1.SetName(treeName)
        myTH1.SetTitle(treeName)

    sigTH1.Sumw2()
    myTH1.Sumw2()

    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetOptStat(0)
        
    rt.gStyle.SetPadLeftMargin(0.15)
    
    c = rt.TCanvas('c','c',500,400)
    sigTH1.Divide(myTH1)
    sigTH1.SetMarkerStyle(8)
    sigTH1.SetMarkerSize(0.8)
    sigTH1.Draw("pe")
    sigTH1.SetMinimum(0)
    sigTH1.GetXaxis().SetTitle('Bin Number')
    sigTH1.GetYaxis().SetTitle('Signal Contamination in %s'%box)
    sigTH1.GetYaxis().SetTitleOffset(2.)

    tleg = rt.TLegend(0.2,0.69,0.5,0.89)
    tleg.AddEntry(sigTH1,'%s (%i, %i)'%(options.model,options.mGluino,options.mLSP))
    tleg.SetLineColor(rt.kWhite)
    tleg.Draw()
    c.Print("%s/signalContamination_%s_%i_%i_%s.pdf"%(options.outDir,options.model,options.mGluino,options.mLSP,box))
    c.Print("%s/signalContamination_%s_%i_%i_%s.C"%(options.outDir,options.model,options.mGluino,options.mLSP,box))
