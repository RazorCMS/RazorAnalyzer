from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import rootTools
from framework import Config


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T2tt",type="string",
                  help="signal model name")
    parser.add_option('--mGluino',dest="mGluino", default=-1,type="float",
                  help="mass of gluino")
    parser.add_option('--mStop',dest="mStop", default=-1,type="float",
                  help="mass of stop")
    parser.add_option('--mLSP',dest="mLSP", default=-1,type="float",
                  help="mass of LSP")
    parser.add_option('-l','--lumi',dest="lumi", default="4000",type="string",
                  help="lumi in pb^-1")
    parser.add_option('-i','--input',dest="inputFile", default="higgsCombineTest.MaxLikelihoodFit.mH120.root",type="string",
                  help="input file with workspace")
    parser.add_option('-f','--fit',dest="fitFile", default="mlfit.root",type="string",
                  help="input file with fit")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store plots")


    (options,args) = parser.parse_args()

    box = options.box
    model = options.model
    
    try: 
        os.environ['CMSSW_BASE']
        loadVal = rt.gSystem.Load("${CMSSW_BASE}/lib/${SCRAM_ARCH}/libHiggsAnalysisCombinedLimit.so")
        if loadVal == -1:
            print "WARNING: NO HIGGS LIBRARY"
    except:
        sys.exit()
    
    if options.mGluino>-1:
        massPoint = 'mGl-%i_mLSP-%i'%(options.mGluino,options.mLSP)
    elif options.mStop>-1:
        massPoint = 'mStop-%i_mLSP-%i'%(options.mStop,options.mLSP)

    
    
    inputFile = rt.TFile(options.inputFile)
    w = inputFile.Get('w')
    w.Print('v')
    fitFile = rt.TFile(options.fitFile)
    fit_b = fitFile.Get('fit_b')
    fit_s = fitFile.Get('fit_s')
    
                                     
    th1x = w.var('th1x')
    fit_b.Print('v')
    fit_s.Print('v')
    
    c = rt.TCanvas('c','c',500,400)
    th1xFrame = th1x.frame()
    
    print "\nBACKGROUND"
    for p in rootTools.RootIterator.RootIterator(fit_b.floatParsFinal()):
        w.var(p.GetName()).setVal(p.getVal())
        w.var(p.GetName()).setError(p.getError())
        print "INITIALIZE PARAMETER %s = %f +- %f"%(p.GetName(),p.getVal(),p.getError())
    w.var('r').setVal(0.)
        
    w.data('data_obs').plotOn(th1xFrame,rt.RooFit.Invisible())
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrame,rt.RooFit.LineColor(rt.kBlack))
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrame,rt.RooFit.Components('shapeBkg_%s_TTj1b_%s'%(box,box)),rt.RooFit.LineColor(rt.kBlue))
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrame,rt.RooFit.Components('shapeBkg_%s_TTj2b_%s'%(box,box)),rt.RooFit.LineColor(rt.kViolet))
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrame,rt.RooFit.Components('shapeBkg_%s_TTj3b_%s'%(box,box)),rt.RooFit.LineColor(rt.kGreen))
    w.data('data_obs').plotOn(th1xFrame)
    th1xFrame.Draw()
    
    c.Print('th1xFrame_b.pdf')

    th1xFrameClone = th1xFrame.emptyClone("th1xFrameClone")
    print "\nSIGNAL"
    for p in rootTools.RootIterator.RootIterator(fit_s.floatParsFinal()):
        w.var(p.GetName()).setVal(p.getVal())
        w.var(p.GetName()).setError(p.getError())
        print "INITIALIZE PARAMETER %s = %f +- %f"%(p.GetName(),p.getVal(),p.getError())
        
    w.data('data_obs').plotOn(th1xFrameClone,rt.RooFit.Invisible())
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrameClone,rt.RooFit.LineColor(rt.kBlack))
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrameClone,rt.RooFit.Components('shapeBkg_%s_TTj1b_%s'%(box,box)),rt.RooFit.LineColor(rt.kBlue))
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrameClone,rt.RooFit.Components('shapeBkg_%s_TTj2b_%s'%(box,box)),rt.RooFit.LineColor(rt.kViolet))
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrameClone,rt.RooFit.Components('shapeBkg_%s_TTj3b_%s'%(box,box)),rt.RooFit.LineColor(rt.kGreen))
    w.pdf('pdf_bin%s_nuis'%box).plotOn(th1xFrameClone,rt.RooFit.Components('shapeSig_%s_%s_%sPdf'%(box,model,box)),rt.RooFit.LineColor(rt.kRed))
    w.data('data_obs').plotOn(th1xFrameClone)
    th1xFrameClone.Draw()
    
    c.Print('th1xFrame_s.pdf')
    
    

