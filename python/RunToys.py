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
    
        
def getTree(myTree,paramNames,nBins,box):
    # first structure
    stringMyStruct1= "{struct MyStruct1{float n2dll_%s;"%box
    for paramName in paramNames:
        stringMyStruct1 = stringMyStruct1+"float %s;" %(paramName)
    for iBinX in range(0,nBins):
        stringMyStruct1 = stringMyStruct1+"float b%i;" %(iBinX)
    #print stringMyStruct1+"};}"
    rando = random.randint(1,999999)
    tempMacro = open("tempMacro_%d.C"%rando,"w")
    tempMacro.write(stringMyStruct1+"};}")
    tempMacro.close()
    rt.gROOT.Macro("tempMacro_%d.C"%rando)
    #rt.gROOT.ProcessLine(stringMyStruct1+"};")
    from ROOT import MyStruct1

    # fills the bins list and associate the bins to the corresponding variables in the structure
    s1 = MyStruct1()
    myTree.Branch('n2dll_%s'%box , rt.AddressOf(s1,'n2dll_%s'%box),'n2dll_%s/F' %box)
    for paramName in paramNames:
        myTree.Branch(paramName , rt.AddressOf(s1,paramName),'%s/F' %paramName)
    for ix in range(0, nBins):
        myTree.Branch("b%i" %(ix) , rt.AddressOf(s1,"b%i" %(ix)),'b%i/F' %ix)

    os.system("rm tempMacro_%d.C"%rando)
    return s1


def runToys(w,options,cfg):
    
    setStyle()
    rt.RooRandom.randomGenerator().SetSeed(options.seed)
    
    extRazorPdf = w.pdf('extRazorPdf')
    dataHist = w.data("data_obs")
    fr = w.obj('fitresult_extRazorPdf_data_obs')
    #nll = w.function('nll_extRazorPdf_data_obs')
    nll = extRazorPdf.createNLL(dataHist)
    th1x = w.var("th1x")
    
    params = extRazorPdf.getParameters(dataHist)
    paramsToRemove = []
    for p in rootTools.RootIterator.RootIterator(params):
        if p.isConstant(): paramsToRemove.append(p)

    [params.remove(p) for p in paramsToRemove]
    paramNames = [p.GetName() for p in rootTools.RootIterator.RootIterator(params)]
    paramNames.sort()
    #params.remove(w.var('BtagCut_TTj1b'))
    #params.remove(w.var('BtagCut_TTj2b'))
    #params.remove(w.var('BtagCut_TTj3b'))
    #params.remove(w.var('MRCut_%s'%options.box))
    #params.remove(w.var('RCut_%s'%options.box))
    
    hessePdf = fr.createHessePdf(params)
    minNll = rt.RooRealVar('minNll_%s'%options.box,'minNll_%s'%options.box,fr.minNll())
    minNll.setConstant(True)
    n2dll = rt.RooFormulaVar('n2dll_%s'%options.box,'2*@0-2*@1',rt.RooArgList(nll,minNll))
    
    #rootTools.Utils.importToWS(w,n2dll)
    #rootTools.Utils.importToWS(w,hessePdf)
    
    x = array('d', cfg.getBinning(options.box)[0]) # MR binning
    y = array('d', cfg.getBinning(options.box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(options.box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    th1x.setBins(nBins)


    unc = 'Bayes'
    if options.freq: unc = 'Freq'
    
    output = rt.TFile.Open(options.outDir+'/toys_%s_%s.root'%(unc,options.box),'recreate')
    output.cd()
    myTree = rt.TTree("myTree", "myTree")
    
    s1 = getTree(myTree,paramNames, nBins,options.box)
    
    value = setattr(s1, 'n2dll_%s'%options.box, n2dll.getVal())
    for p in rootTools.RootIterator.RootIterator(fr.floatParsFinal()):
        w.var(p.GetName()).setVal(p.getVal())
        w.var(p.GetName()).setError(p.getError())
        value = setattr(s1, p.GetName(), p.getVal())
        
    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())
    
    for iBinX in range(0,nBins):
        th1x.setVal(iBinX+0.5)
        value = setattr(s1, 'b%i'%iBinX, float(asimov.weight(rt.RooArgSet(th1x))))
    
    myTree.Fill()
        
    iToy = 0
    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)
    
    nBadPars = 0
    while iToy < options.nToys:
        if options.freq:
            pSet = fr.floatParsFinal()
        else:
            if options.noSys:                
                pSet = fr.floatParsFinal()
            else:                
                pSet = fr.randomizePars()
        for p in rootTools.RootIterator.RootIterator(pSet):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())

        badPars = []
        for bkgd in ['TTj1b','TTj2b']:
            badPars.append(w.var('n_%s_%s'%(bkgd,options.box)).getVal() <= 0)
            badPars.append(w.var('b_%s_%s'%(bkgd,options.box)).getVal() <= 0)
            badPars.append(w.var('MR0_%s_%s'%(bkgd,options.box)).getVal() >= w.var('MR').getMin())
            badPars.append(w.var('R0_%s_%s'%(bkgd,options.box)).getVal()  >=  w.var('Rsq').getMin())
        if any(badPars):
            nBadPars+=1
            #print "BAD PARS #%i"%nBadPars
            continue
        
        errorCountBefore = rt.RooMsgService.instance().errorCount()
        pdfValV = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
        errorCountAfter = rt.RooMsgService.instance().errorCount()
        if errorCountAfter > errorCountBefore: continue

        errorCountBefore = rt.RooMsgService.instance().errorCount()
        n2dll.getVal()
        errorCountAfter = rt.RooMsgService.instance().errorCount()
        if errorCountAfter > errorCountBefore: continue

        #if sigmaCut>0 and math.sqrt(n2dll.getVal())>sigmaCut: continue
        
        errorCountBefore = rt.RooMsgService.instance().errorCount()
        if options.noStat:
            asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy_%i'%iToy),rt.RooFit.Asimov())
        else:
            asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy_%i'%iToy))
        errorCountAfter = rt.RooMsgService.instance().errorCount()
        if errorCountAfter > errorCountBefore: continue

        pSetSave = pSet
        if options.freq:
            fr_toy = extRazorPdf.fitTo(asimov,rt.RooFit.Save(),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.Warnings(False))
            if fr_toy.covQual() < 3: continue
            pSetSave = fr_toy.floatParsFinal()
            if options.noStat:
                asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('pred_%i'%iToy),rt.RooFit.Asimov())
            else:
                asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('pred_%i'%iToy))

        value = setattr(s1, 'n2dll_%s'%options.box, n2dll.getVal())

        for p in rootTools.RootIterator.RootIterator(pSetSave):
            value = setattr(s1, p.GetName(), p.getVal())

        for iBinX in range(0,nBins):
            th1x.setVal(iBinX+0.5)
            pdfValV = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
            #if iToy%1000==0: print iBinX, asimov.weight(rt.RooArgSet(th1x)), pdfValV
            value = setattr(s1, 'b%i'%iBinX, float(asimov.weight(rt.RooArgSet(th1x))))
        
        if iToy%100==0: print "generating %ith toy"%iToy
        myTree.Fill()
        iToy+=1
    rt.RooMsgService.instance().reset()

    
    w.Print('v')
    output.cd()
    myTree.Write()
    w.Write()
    output.Close()
    return output.GetName()

        
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
    parser.add_option('-t','--toys',dest="nToys", default=3000,type="int",
                  help="number of toys")
    parser.add_option('-s','--seed',dest="seed", default=1988,type="int",
                  help="random seed")
    parser.add_option('--no-stat',dest="noStat",default=False,action='store_true',
                  help="no statistical uncertainty, just systematic uncertainty, default is statstical + systematic uncertainty")
    parser.add_option('--no-sys',dest="noSys",default=False,action='store_true',
                  help="no systematic uncertainty, just statistical uncertainty, default is statstical + systematic uncertainty")
    parser.add_option('--freq',dest="freq",default=False,action='store_true',
                  help="refit each toy with only statistical fluctuations, as in frequentist approach; default is bayeseian")
    parser.add_option('-i','--inputFitFile',dest="inputFitFile", default=None,type="string",
                  help="input fit file")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    
    lumi = options.lumi

    inputFitFile = options.inputFitFile
    
    
    lumi_in = 0.

    if inputFitFile is not None:
        rootFile = rt.TFile.Open(inputFitFile,"r")
        w = rootFile.Get("w"+options.box)
        fr = w.obj("fitresult_extRazorPdf_data_obs")
        nll = w.function("nll_extRazorPdf_data_obs")
    else:
        for f in args:
            if f.lower().endswith('.root'):
                rootFile = rt.TFile(f)
                workspace = rootFile.Get('w'+options.box)
                data = workspace.data('RMRTree')
                lumi_in = 1000.*float([g.replace('lumi-','') for g in f.split('_') if g.find('lumi')!=-1][0])

  
        w = rt.RooWorkspace("w"+options.box)
    
        paramNames, bkgs = initializeWorkspace(w,cfg,options.box)
        rootTools.Utils.importToWS(w,data)
    
        th1x = w.var('th1x')
    
        myTH1 = convertDataset2TH1(data, cfg, options.box, w)
        myTH1.Scale(lumi/lumi_in)
        dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), myTH1)
        rootTools.Utils.importToWS(w,dataHist)

        nll = extRazorPdf.createNLL(dataHist)
        #m = rt.RooMinuit(nll)
        #m.migrad()
        #m.hesse()
        #fr = m.save()
        fr = extRazorPdf.fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Minimizer('Minuit2','improve'),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.SumW2Error(False))
        fr.Print('v')
        rootTools.Utils.importToWS(w,fr)
        #rootTools.Utils.importToWS(w,nll)

        
    outputName = runToys(w,options,cfg)

    print "writing tree to %s"%(outputName)
