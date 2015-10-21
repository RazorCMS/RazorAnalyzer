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
    
    rando = random.randint(1,999999)
    # first structure
    stringMyStruct1= "void tempMacro_%d(){struct MyStruct1{float n2llr_%s; float chi2_%s; int covQual_%s;"%(rando,box,box,box)
    for paramName in paramNames:
        stringMyStruct1 = stringMyStruct1+"float %s; float %s_error;" %(paramName,paramName)
    for iBinX in range(0,nBins):
        stringMyStruct1 = stringMyStruct1+"float b%i;" %(iBinX)
    #print stringMyStruct1+"};}"
    tempMacro = open("tempMacro_%d.C"%rando,"w")
    tempMacro.write(stringMyStruct1+"};}")
    tempMacro.close()
    rt.gROOT.Macro("tempMacro_%d.C"%rando)
    #rt.gROOT.ProcessLine(".x tempMacro_%d.C"%rando)
    #rt.gROOT.ProcessLine(stringMyStruct1+"};")
    from ROOT import MyStruct1

    # fills the bins list and associate the bins to the corresponding variables in the structure
    s1 = MyStruct1()
    myTree.Branch('n2llr_%s'%box , rt.AddressOf(s1,'n2llr_%s'%box),'n2llr_%s/F' %box)
    myTree.Branch('chi2_%s'%box , rt.AddressOf(s1,'chi2_%s'%box),'chi2_%s/F' %box)
    myTree.Branch('covQual_%s'%box , rt.AddressOf(s1,'covQual_%s'%box),'covQual_%s/I' %box)
    for paramName in paramNames:
        myTree.Branch(paramName , rt.AddressOf(s1,paramName),'%s/F' %paramName)
        myTree.Branch('%s_error'%paramName , rt.AddressOf(s1,'%s_error'%paramName),'%s_error/F' %paramName)
    for ix in range(0, nBins):
        myTree.Branch("b%i" %(ix) , rt.AddressOf(s1,"b%i" %(ix)),'b%i/F' %ix)

    #os.system("rm tempMacro_%d.C"%rando)
    return s1


def runToys(w,options,cfg):
    
    setStyle()
    rt.RooRandom.randomGenerator().SetSeed(options.seed)
    
    extRazorPdf = w.pdf('extRazorPdf')
    dataHist = w.data("data_obs")        
    fr = w.obj('fitresult_extRazorPdf_data_obs')

    if options.r>-1:
        extSpBPdf = w.pdf('extSpBPdf')
                
    #nll = w.function('nll_extRazorPdf_data_obs')
    nll = extRazorPdf.createNLL(dataHist,rt.RooFit.Extended(True))
    chi2 = extRazorPdf.createChi2(dataHist,rt.RooFit.Extended(True))
    th1x = w.var("th1x")
    
    params = extRazorPdf.getParameters(dataHist)
    paramsToRemove = []
    for p in rootTools.RootIterator.RootIterator(params):
        if p.isConstant(): paramsToRemove.append(p)

    [params.remove(p) for p in paramsToRemove]
    paramNames = [p.GetName() for p in rootTools.RootIterator.RootIterator(params)]
    paramNames.sort()    
    if options.r>-1: paramNames.append('r')
    #params.remove(w.var('BtagCut_TTj1b'))
    #params.remove(w.var('BtagCut_TTj2b'))
    #params.remove(w.var('BtagCut_TTj3b'))
    #params.remove(w.var('MRCut_%s'%options.box))
    #params.remove(w.var('RCut_%s'%options.box))
    
    #hessePdf = fr.createHessePdf(params)
    #minNll = rt.RooRealVar('minNll_%s'%options.box,'minNll_%s'%options.box,fr.minNll())
    #minNll.setConstant(True)
    #n2llr = rt.RooFormulaVar('n2llr_%s'%options.box,'2*@0-2*@1',rt.RooArgList(nll,minNll))
    
    #rootTools.Utils.importToWS(w,n2llr)
    #rootTools.Utils.importToWS(w,hessePdf)
    
    x = array('d', cfg.getBinning(options.box)[0]) # MR binning
    y = array('d', cfg.getBinning(options.box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(options.box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    th1x.setBins(nBins)


    unc = 'Bayes'
    if options.noStat: unc = "Bayes_noStat"
    if options.noSys: unc = "Bayes_noSys"
    if options.freq: unc = 'Freq'

    if options.r>-1:
        rString = str('%.3f'%options.r).replace(".","p")
        output = rt.TFile.Open(options.outDir+'/toys_%s_r%s_%s.root'%(unc,rString,options.box),'recreate')
    else:
        output = rt.TFile.Open(options.outDir+'/toys_%s_%s.root'%(unc,options.box),'recreate')
        
    output.cd()
    myTree = rt.TTree("myTree", "myTree")
    
    s1 = getTree(myTree,paramNames, nBins,options.box)
    
    for p in rootTools.RootIterator.RootIterator(fr.floatParsFinal()):
        w.var(p.GetName()).setVal(p.getVal())
        w.var(p.GetName()).setError(p.getError())
        value = setattr(s1, p.GetName(), p.getVal())
        value = setattr(s1, p.GetName()+'_error', p.getError())
    if options.r>-1:
        value =  setattr(s1, 'r', 0)
        value =  setattr(s1, 'r_error', 0)
        
    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('central'),rt.RooFit.Asimov())

    chi2_data = 0
    n2llr_data = 0
    for iBinX in range(0,nBins):
        th1x.setVal(iBinX+0.5)
        expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
        observed = float(dataHist.weight(rt.RooArgSet(th1x)))
        value = setattr(s1, 'b%i'%iBinX, expected)
        if expected>0:
            chi2_data += ( observed - expected ) * ( observed - expected ) / ( expected )
        if observed>0:
            n2llr_data += 2 * ( observed*rt.TMath.Log(observed/expected) - observed )
        n2llr_data += 2 * ( expected )

        
    value = setattr(s1, 'n2llr_%s'%options.box, n2llr_data)
    value = setattr(s1, 'chi2_%s'%options.box, chi2_data)
    myTree.Fill()
        
    iToy = 0
    
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
        for bkgd in ['TTj0b','TTj1b','TTj2b','TTj3b']:
            if w.var('n_%s_%s'%(bkgd,options.box))!=None:                
                badPars.append(w.var('n_%s_%s'%(bkgd,options.box)).getVal() <= 0)
            if w.var('b_%s_%s'%(bkgd,options.box))!=None:                                
                badPars.append(w.var('b_%s_%s'%(bkgd,options.box)).getVal() <= 0)
            if w.var('MR0_%s_%s'%(bkgd,options.box))!=None:                                
                badPars.append(w.var('MR0_%s_%s'%(bkgd,options.box)).getVal() >= w.var('MR').getMin())
            if w.var('R0_%s_%s'%(bkgd,options.box))!=None:                   
                badPars.append(w.var('R0_%s_%s'%(bkgd,options.box)).getVal()  >=  w.var('Rsq').getMin())
        if any(badPars):
            nBadPars+=1
            #print "bad pars toy=%i"%iToy
            continue

        #print "good pars"        
        #for bkgd in ['TTj0b','TTj1b','TTj2b','TTj3b']:
        #    if w.var('n_%s_%s'%(bkgd,options.box))!=None:                
        #        print 'n_%s_%s = %f'%(bkgd,options.box,w.var('n_%s_%s'%(bkgd,options.box)).getVal())
        #    if w.var('b_%s_%s'%(bkgd,options.box))!=None:                                
        #        print 'b_%s_%s = %f'%(bkgd,options.box,w.var('b_%s_%s'%(bkgd,options.box)).getVal())
        #    if w.var('MR0_%s_%s'%(bkgd,options.box))!=None:                                
        #        print 'MR0_%s_%s = %f'%(bkgd,options.box,w.var('MR0_%s_%s'%(bkgd,options.box)).getVal())
        #    if w.var('R0_%s_%s'%(bkgd,options.box))!=None:                   
        #        print 'R0_%s_%s = %f'%(bkgd,options.box,w.var('R0_%s_%s'%(bkgd,options.box)).getVal())
                
        errorCountBefore = rt.RooMsgService.instance().errorCount()
        th1x.setVal(0.5) # check number of events in first bin
        pdfValV = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
        pdfVal0 = extRazorPdf.getValV(0) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
        errorCountAfter = rt.RooMsgService.instance().errorCount()
        if errorCountAfter > errorCountBefore:            
            #print "can't evaulate pdf toy=%i"%iToy
            continue
        
        if pdfValV<=0.0 and pdfVal0<=0.0:
            #print "pdf valv = %f"%(pdfValV)
            #print "pdf val0 = %f"%(pdfVal0)
            continue
        
        errorCountBefore = rt.RooMsgService.instance().errorCount()        
        #print "start generating toy=%i"%iToy
        if options.noStat:
            #asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy_%i'%iToy),rt.RooFit.Asimov())            
            if options.r>-1:
                w.var('r').setVal(options.r)            
                asimov = extSpBPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Asimov())
            else:
                asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'),rt.RooFit.Asimov())
        else:
            #asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy_%i'%iToy))
            if options.r>-1:                
                w.var('r').setVal(options.r)      
                asimov = extSpBPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'))
            else:
                asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy'))

        #print "toy entries = %i"%asimov.sumEntries()
        errorCountAfter = rt.RooMsgService.instance().errorCount()   
        if errorCountAfter > errorCountBefore:
            #print "can't generate toy=%i"%iToy
            continue

        #print "SUCCESS: generated toy=%i"%iToy

        pSetSave = pSet
        if options.freq:                      
            if options.r>-1:
                fr_toy = extSpBPdf.fitTo(asimov,rt.RooFit.Save(),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.Warnings(False))
            else:
                fr_toy = extRazorPdf.fitTo(asimov,rt.RooFit.Save(),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.Warnings(False))
            #if fr_toy.covQual() < 2: continue            
            value = setattr(s1,'covQual_%s'%options.box, fr_toy.covQual())
            pSetSave = fr_toy.floatParsFinal()
            #if options.noStat:
            #    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('pred_%i'%iToy),rt.RooFit.Asimov())
            #else:
            #    asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('pred_%i'%iToy))

                
        for p in rootTools.RootIterator.RootIterator(pSetSave):
            value = setattr(s1, p.GetName(), p.getVal())
            value = setattr(s1, p.GetName()+"_error", p.getError())

        chi2_toy = 0
        n2llr_toy = 0
        # restore best-fit to calculate expected values        
        for p in rootTools.RootIterator.RootIterator(pSetSave):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())
        for iBinX in range(0,nBins):
            th1x.setVal(iBinX+0.5)
            expected = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
            toy = float(asimov.weight(rt.RooArgSet(th1x)))
            observed = float(dataHist.weight(rt.RooArgSet(th1x)))
            value = setattr(s1, 'b%i'%iBinX, toy)
            if expected>0:
                chi2_toy += ( toy - expected ) * ( toy - expected ) / ( expected )
            if toy>0 and expected>0:
                n2llr_toy += 2 * ( toy*rt.TMath.Log(toy/expected) - toy )
            n2llr_toy += 2 * ( expected )
            
        value = setattr(s1, 'n2llr_%s'%options.box, n2llr_toy)
        value = setattr(s1, 'chi2_%s'%options.box, chi2_toy)
        
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
    parser.add_option('-r','--signal-strength',dest="r", default=-1,type="float",
                  help="signal strength => do each fit the the SpB pdf")
    parser.add_option('-i','--input-fit-file',dest="inputFitFile", default=None,type="string",
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

        #m = rt.RooMinuit(nll)
        #m.migrad()
        #m.hesse()
        #fr = m.save()
        fr = extRazorPdf.fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Minimizer('Minuit2','improve'),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.SumW2Error(False))
        fr.Print('v')
        rootTools.Utils.importToWS(w,fr)

        
    outputName = runToys(w,options,cfg)

    print "writing tree to %s"%(outputName)
