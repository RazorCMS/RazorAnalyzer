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
    
def find68ProbRange(hToy, probVal=0.6827):
    minVal = 0.
    maxVal = 100000.
    if hToy.Integral()<=0: return hToy.GetBinCenter(hToy.GetMaximumBin()),max(minVal,0.),maxVal
    # get the bin contents
    probsList = []
    for  i in range(1, hToy.GetNbinsX()+1):
        probsList.append(hToy.GetBinContent(i)/hToy.Integral())
    probsList.sort()
    prob = 0
    prob68 = 0
    found = False
    for i in range(0,len(probsList)):
        if prob+probsList[i] >= 1-probVal and not found:
            frac = (1-probVal-prob)/probsList[i]
            prob68 = probsList[i-1]+frac*(probsList[i]-probsList[i-1])
            found = True
        prob = prob + probsList[i]

    foundMin = False
    foundMax = False
    for  i in range(0, hToy.GetNbinsX()):
        if not foundMin and hToy.GetBinContent(i+1) >= prob68:
            fraction = (prob68-hToy.GetBinContent(i))/(hToy.GetBinContent(i+1)-hToy.GetBinContent(i))
            minVal = hToy.GetBinLowEdge(i)+hToy.GetBinWidth(i)*fraction
            foundMin = True
        if not foundMax and hToy.GetBinContent(hToy.GetNbinsX()-i) >= prob68:
            fraction = (prob68-hToy.GetBinContent(hToy.GetNbinsX()-i+1))/(hToy.GetBinContent(hToy.GetNbinsX()-i)-hToy.GetBinContent(hToy.GetNbinsX()-i+1))
            maxVal = hToy.GetBinLowEdge(hToy.GetNbinsX()-i)+hToy.GetBinWidth(hToy.GetNbinsX()-i)*(1-fraction)
            foundMax = True
    return hToy.GetBinCenter(hToy.GetMaximumBin()),max(minVal,0.),maxVal

def findMedian(myHisto):
    prob = 0
    median = 0
    for i in range(1, myHisto.GetNbinsX()+1):
        if prob <= 0.5 and prob+myHisto.GetBinContent(i) > 0.5:
            median = myHisto.GetBinCenter(i)
        prob = prob + myHisto.GetBinContent(i)
    return median

def getPValue(n, hToy):
    if hToy.Integral() <= 0.: return 0.
    Prob_n = hToy.GetBinContent(hToy.FindBin(n))
    Prob = 0
    for i in range(1, hToy.GetNbinsX()+1):
        if hToy.GetBinContent(i)<= Prob_n: Prob += hToy.GetBinContent(i)
    Prob = Prob/hToy.Integral()
    return Prob

def getSigmaFromPval(n, hToy, pVal):
    if hToy.GetMaximumBin() == hToy.FindBin(n): return 0.
    medianVal = findMedian(hToy)
    coreProb = 1. - pVal
    if n>medianVal: return rt.TMath.NormQuantile(0.5 + coreProb/2.)
    else: return -rt.TMath.NormQuantile(0.5 + coreProb/2.)
        
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

def getBinSumDicts(sumType, minX, maxX, minY, maxY, minZ, maxZ, x, y, z):

    binSumDict = {}

    if sumType.lower()=="x":
        for i in range(1,len(x)):
            binXSumDict[i] = ''
    elif sumType.lower()=="y":
        for j in range(1,len(y)):
            binYSumDict[j] = ''
    elif sumType.lower()=="yx":
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                binSumDict[i,j] = ''
    elif sumType.lower()=="zyx":
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                for k in range(1,len(z)):
                    binSumDict[i,j,k] = ''
            
    iBinX = -1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                if i>=minX and i<=maxX and j>=minY and j<=maxY and z>=minZ and i<=maxZ:
                    if sumType.lower()=="x":
                        binSumDict[i] += 'b%i+'%iBinX
                    elif sumType.lower()=="y":
                        binSumDict[j] += 'b%i+'%iBinX
                    elif sumType.lower()=="yx":
                        binSumDict[i,j] += 'b%i+'%iBinX
                    elif sumType.lower()=="zyx":
                        binSumDict[i,j,k] = 'b%i'%iBinX
                
    if sumType.lower()=="x":
        for i in range(1,len(x)):
            binSumDict[i] = binSumDict[i][0:-1]            
    elif sumType.lower()=="y":
        for j in range(1,len(y)):
            binSumDict[j] = binSumDict[j][0:-1]      
    elif sumType.lower()=="yx":
        for i in range(1,len(x)):
            for j in range(1,len(y)):
                binSumDict[i,j] = binSumDict[i,j][0:-1]
            
    return binSumDict

def getBestFitRms(myTree, sumName, nObs, c, outDir):
    myTree.GetEntry(0)
    bestFit = eval(sumName.replace('b','myTree.b'))
    myTree.Draw('%s>>htest%s'%(sumName,sumName.replace('+','')))
    htemp = rt.gPad.GetPrimitive("htest%s"%sumName.replace('+',''))
    #xmax = int(htemp.GetXaxis().GetXmax()+1)
    xmax = max(htemp.GetXaxis().GetXmax(),nObs)
    xmin = max(0,htemp.GetXaxis().GetXmin())
    nbins = htemp.GetNbinsX()
    nbins = 100.
    #xmode, xminus, xplus = find68ProbRange(htemp, probVal=0.5)
    #binWidth = 2.*(xplus-xminus)*rt.TMath.Power(htemp.GetEntries(),-1./3)
    #binWidth = max(1.,binWidth)
    #nbins = int((xmax-xmin)/binWidth)
    myTree.Draw('%s>>htemp%s(%i,%f,%f)'%(sumName,sumName.replace('+',''),nbins,xmin,xmax))
    htemp = rt.gPad.GetPrimitive("htemp%s"%sumName.replace('+',''))
    mean = htemp.GetMean()
    mode = htemp.GetBinCenter(htemp.GetMaximumBin())
    rms = htemp.GetRMS()
    pvalue = getPValue(nObs,htemp)
    nsigma = getSigmaFromPval(nObs,htemp,pvalue)
    print '%s, bestFit %f, mean %.1f, mode %.1f, rms %.1f, pvalue %f, nsigma %.1f'%(sumName, bestFit,mean,mode,rms,pvalue,nsigma)

    htemp.SetMaximum(1.1*(htemp.GetMaximum()+htemp.GetBinError(htemp.GetMaximumBin())))
    htemp.SetMarkerStyle(20)
    htemp.SetMarkerColor(rt.kBlack)
    htemp.SetLineColor(rt.kBlack)
    tgraph = rt.TGraph(4)
    tgraph.SetPoint(1, bestFit-rms,0)
    tgraph.SetPoint(2, bestFit-rms,htemp.GetMaximum())
    tgraph.SetPoint(3, bestFit+rms,htemp.GetMaximum())
    tgraph.SetPoint(4, bestFit+rms,0)
    tgraph.SetFillColor(rt.kBlue-10)
    tgraph.SetLineColor(rt.kBlue-10)
    tline = rt.TLine(bestFit,0,bestFit,htemp.GetMaximum())
    tline.SetLineColor(rt.kBlue)
    tline.SetLineWidth(2)
    tlineObs = rt.TLine(nObs,0,nObs,htemp.GetMaximum())
    tlineObs.SetLineColor(rt.kBlack)
    tlineObs.SetLineWidth(2)
    htemp.Draw("pe")
    tgraph.Draw("fsame")
    tline.Draw("same")
    tlineObs.Draw("same")
    htemp.Draw("pesame")
    htemp.Fit("gaus")
    c.Print(outDir+'/%s.pdf'%sumName.replace('+',''))
    c.Print(outDir+'/%s.C'%sumName.replace('+',''))
    return bestFit, rms, pvalue, nsigma


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
    sys.exit()
    
    output = rt.TFile.Open(outputName,"read")
    myTree = output.Get("myTree")


    
    c = rt.TCanvas('c','c',500,400)
    
    x = array('d', cfg.getBinning(options.box)[0]) # MR binning
    y = array('d', cfg.getBinning(options.box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(options.box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    rt.TH1D.SetDefaultSumw2()
    rt.TH2D.SetDefaultSumw2()
    rt.TH3D.SetDefaultSumw2()
    h_nBtagRsqMR = rt.TH3D("h_nBtagRsqMR","h_nBtagRsqMR",len(x)-1,x,len(y)-1,y,len(z)-1,z)
    h_RsqMR_ = {}; h_MR_ = {}; h_Rsq_ = {}; h_FF_ = {};
    for btag in [[1,2,3],[1],[2],[3]]:
        btagString = [str(k) for k in btag]
        btagString = "_".join(btagString)
        print btagString
        h_RsqMR_[btagString] = rt.TH2D("h_RsqMR_%s"%btagString,"h_RsqMR_%s"%btagString,len(x)-1,x,len(y)-1,y)
        h_FF_[btagString] = rt.TH2D("h_FF_%s"%btagString,"h_FF_%s"%btagString,len(x)-1,x,len(y)-1,y)
        h_MR_[btagString] = rt.TH1D("h_MR_%s"%btagString,"h_MR_%s"%btagString,len(x)-1,x)
        h_Rsq_[btagString] = rt.TH1D("h_Rsq_%s"%btagString,"h_Rsq_%s"%btagString,len(y)-1,y)

    h_data_th1x = dataHist.createHistogram('h_data_th1x',th1x)
    h_data_nBtagRsqMR = rt.TH3D("h_data_nBtagRsqMR","h_data_nBtagRsqMR",len(x)-1,x,len(y)-1,y,len(z)-1,z)
    h_data_RsqMR_ = {}; h_data_MR_ = {}; h_data_Rsq_ = {}
    
    iBinX = -1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                h_data_nBtagRsqMR.SetBinContent(i,j,k,h_data_th1x.GetBinContent(iBinX+1))
                h_data_nBtagRsqMR.SetBinError(i,j,k,h_data_th1x.GetBinError(iBinX+1))
    h_data_nBtagRsqMR_ = {}
    h_data_nBtagRsqMR_['1'] =  h_data_nBtagRsqMR.Clone(h_data_nBtagRsqMR.GetName()+'_1')
    h_data_nBtagRsqMR_['2'] =  h_data_nBtagRsqMR.Clone(h_data_nBtagRsqMR.GetName()+'_2')
    h_data_nBtagRsqMR_['3'] =  h_data_nBtagRsqMR.Clone(h_data_nBtagRsqMR.GetName()+'_3')
                                            
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            h_data_nBtagRsqMR_['1'].SetBinContent(i,j,2,0)
            h_data_nBtagRsqMR_['1'].SetBinContent(i,j,3,0)
            h_data_nBtagRsqMR_['2'].SetBinContent(i,j,1,0)
            h_data_nBtagRsqMR_['2'].SetBinContent(i,j,3,0)
            h_data_nBtagRsqMR_['3'].SetBinContent(i,j,1,0)
            h_data_nBtagRsqMR_['3'].SetBinContent(i,j,2,0)
            h_data_nBtagRsqMR_['1'].SetBinError(i,j,2,0)
            h_data_nBtagRsqMR_['1'].SetBinError(i,j,3,0)
            h_data_nBtagRsqMR_['2'].SetBinError(i,j,1,0)
            h_data_nBtagRsqMR_['2'].SetBinError(i,j,3,0)
            h_data_nBtagRsqMR_['3'].SetBinError(i,j,1,0)
            h_data_nBtagRsqMR_['3'].SetBinError(i,j,2,0)
            

    h_data_RsqMR_['1_2_3'] = h_data_nBtagRsqMR.Project3D("yxe")
    h_data_MR_['1_2_3'] = h_data_nBtagRsqMR.Project3D("xe")
    h_data_Rsq_['1_2_3'] = h_data_nBtagRsqMR.Project3D("ye")

    for btag in ['1','2','3']:
        h_data_RsqMR_[btag] = h_data_nBtagRsqMR_[btag].Project3D("yxe")
        h_data_MR_[btag] = h_data_nBtagRsqMR_[btag].Project3D("xe")
        h_data_Rsq_[btag] = h_data_nBtagRsqMR_[btag].Project3D("ye")

        
    for btag in [[1,2,3],[1],[2],[3]]:
        btagString = [str(k) for k in btag]
        btagString = "_".join(btagString)
        binXSumDict = getBinSumDicts("x",0,len(x),0,len(y),0,len(z),x,y,z)
        binXSumDict, binYSumDict, binXYSumDict, binXYZSumDict = getBinSumDicts("x",0,len(x),0,len(y),0,len(z),x,y,z)
        for i, sumName in binXSumDict.iteritems():
            nObs = h_data_MR_[btagString].GetBinContent(i)
            bestFit, rms, pvalue, nsigma = getBestFitRms(myTree,sumName,nObs,c,options.outDir)
            h_MR_[btagString].SetBinContent(i,bestFit)
            h_MR_[btagString].SetBinError(i,rms)
        for j, sumName in binYSumDict.iteritems():
            nObs = h_data_Rsq_[btagString].GetBinContent(j)
            bestFit, rms, pvalue, nsigma = getBestFitRms(myTree,sumName,nObs,c,options.outDir)
            h_Rsq_[btagString].SetBinContent(j,bestFit)
            h_Rsq_[btagString].SetBinError(j,rms)
        for (i,j), sumName in binXYSumDict.iteritems():
            nObs = h_data_RsqMR_[btagString].GetBinContent(i,j)
            bestFit, rms, pvalue, nsigma = getBestFitRms(myTree,sumName,nObs,c,options.outDir)
            h_RsqMR_[btagString].SetBinContent(i,j,bestFit)
            h_RsqMR_[btagString].SetBinError(i,j,rms)
            if nsigma<0: h_FF_[btagString].SetBinContent(i,j,nsigma)
            
    binXSumDict, binYSumDict, binXYSumDict, binXYZSumDict = getBinSumDicts(x,y,z,btag=[1,2,3])
    for (i,j,k), sumName in binXYZSumDict.iteritems():
        nObs = h_data_nBtagRsqMR.GetBinContent(i,j,k)
        bestFit, rms, pvalue, nsigma = getBestFitRms(myTree,sumName,nObs,c,options.outDir)
        h_nBtagRsqMR.SetBinContent(i,j,k,bestFit)
        h_nBtagRsqMR.SetBinError(i,j,k,rms)
            
    
    unc = 'Bayes'
    if options.freq: unc = 'Freq'
