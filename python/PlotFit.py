from optparse import OptionParser
import os
import ROOT as rt
from array import *
import sys
import rootTools
from framework import Config



def initializeWorkspace(w,cfg,box):
    parameters = cfg.getVariables(box, "combine_parameters")
    paramNames = []
    for parameter in parameters:
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if paramName.find("Cut")==-1:
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
        else:
            w.var(paramName).setConstant(True)
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    w.factory('th1x[0,0,%i]'%nBins)
    emptyHist3D = rt.TH3D("emptyHist3D","emptyHist3D",len(x)-1,x,len(y)-1,y,len(z)-1,z)
    #rootTools.Utils.importToWS(w,emptyHist3D)
    #combine = cfg.getPdfs(box, "combine_pdfs",w)

    w.Print('v')
    razorPdf_TTj1b = rt.RooRazor3DBinPdf("%s_%s"%(box,"TTj1b"),"razorPdf_%s_%s"%(box,"TTj1b"),
                                             w.var("th1x"),
                                             w.var("MR0_%s_%s"%("TTj1b",box)),w.var("R0_%s_%s"%("TTj1b",box)),
                                             w.var("b_%s_%s"%("TTj1b",box)),w.var("n_%s_%s"%("TTj1b",box)),
                                             w.var("MRCut_%s"%box),w.var("RCut_%s"%box),w.var("BtagCut_%s"%("TTj1b")),
                                             emptyHist3D)
    razorPdf_TTj2b = rt.RooRazor3DBinPdf("%s_%s"%(box,"TTj2b"),"razorPdf_%s_%s"%(box,"TTj2b"),
                                             w.var("th1x"),
                                             w.var("MR0_%s_%s"%("TTj2b",box)),w.var("R0_%s_%s"%("TTj2b",box)),
                                             w.var("b_%s_%s"%("TTj2b",box)),w.var("n_%s_%s"%("TTj2b",box)),
                                             w.var("MRCut_%s"%box),w.var("RCut_%s"%box),w.var("BtagCut_%s"%("TTj2b")),
                                             emptyHist3D)
    razorPdf_TTj3b = rt.RooRazor3DBinPdf("%s_%s"%(box,"TTj3b"),"razorPdf_%s_%s"%(box,"TTj3b"),
                                             w.var("th1x"),
                                             w.var("MR0_%s_%s"%("TTj2b",box)),w.var("R0_%s_%s"%("TTj2b",box)),
                                             w.var("b_%s_%s"%("TTj2b",box)),w.var("n_%s_%s"%("TTj2b",box)),
                                             w.var("MRCut_%s"%box),w.var("RCut_%s"%box),w.var("BtagCut_%s"%("TTj3b")),
                                             emptyHist3D)
    rootTools.Utils.importToWS(w,razorPdf_TTj1b)
    rootTools.Utils.importToWS(w,razorPdf_TTj2b)
    rootTools.Utils.importToWS(w,razorPdf_TTj3b)

    return paramNames




if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use"),
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="T1bbbb",type="string",
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
    
    cfg = Config.Config(options.config)
    
    model = options.model
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
        

    th1x.Print()
