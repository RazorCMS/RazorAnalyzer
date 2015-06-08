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

seed = 1988

def setStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)

    
def print1DCanvas(c,h,h_data,printName,xTitle,yTitle,lumiLabel="",boxLabel="",tLeg=None):
    c.SetLogx(0)
    c.SetLogy(1)
    
    pad1 = rt.TPad(c.GetName()+"_pad1","pad1",0,0.25,1,1)
    pad2 = rt.TPad(c.GetName()+"_pad2","pad2",0,0,1,0.25)
    pad1.Range(-213.4588,-0.3237935,4222.803,5.412602);
    pad2.Range(-213.4588,-2.206896,4222.803,3.241379);
    pad1.SetLeftMargin(0.15)
    pad2.SetLeftMargin(0.15)
    pad1.SetRightMargin(0.05)
    pad2.SetRightMargin(0.05)
    pad1.SetTopMargin(0.12)
    pad2.SetTopMargin(0.)
    pad1.SetBottomMargin(0.)
    pad2.SetBottomMargin(0.47)
    pad1.Draw()
    pad1.cd()
    rt.gPad.SetLogy()
    
    h.SetLineWidth(2)
    h.SetLineColor(rt.kBlue)
    hClone = h.Clone(h.GetName()+"Clone")
    hClone.SetLineColor(rt.kBlue)
    hClone.SetFillColor(rt.kBlue-10)
    h_data.SetMarkerColor(rt.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetLineColor(rt.kBlack)
    h_data.GetXaxis().SetTitle(xTitle)
    h_data.GetYaxis().SetTitle(yTitle)
    h_data.GetXaxis().SetLabelOffset(0.16)
    h_data.GetXaxis().SetLabelSize(0.06)
    h_data.GetYaxis().SetLabelSize(0.06)
    h_data.GetXaxis().SetTitleSize(0.06)
    h_data.GetYaxis().SetTitleSize(0.08)
    h_data.GetXaxis().SetTitleOffset(0.8)
    h_data.GetYaxis().SetTitleOffset(0.7)
    h_data.GetXaxis().SetTicks("+-")
    h_data.SetMaximum(math.pow(h_data.GetBinContent(h_data.GetMaximumBin()),1.25))
    h_data.SetMinimum(1e-1*h_data.GetBinContent(h_data.GetMinimumBin()))
    h_data.Draw("pe")
    hClone.Draw("e2same")
    h.SetFillStyle(0)
    h.DrawCopy("histsame")
    h_data.Draw("pesame")
    pad1.Draw()
    c.Update()
    c.cd()
    pad2.Draw()
    pad2.cd()
    rt.gPad.SetLogy(0)
    
    h_data.Sumw2()
    hClone.Sumw2()
    hDataDivide = h_data.Clone(h_data.GetName()+"Divide")
    hDataDivide.Sumw2()

    hDivide = h.Clone(h.GetName()+"Divide") 
    hCloneDivide = hClone.Clone(hClone.GetName()+"Divide") 
    hCloneDivide.GetYaxis().SetLabelSize(0.18)
    hCloneDivide.SetTitle("")
    hCloneDivide.SetMaximum(3.5)
    hCloneDivide.SetMinimum(0.)
    hCloneDivide.GetXaxis().SetLabelSize(0.22)
    hCloneDivide.GetXaxis().SetTitleSize(0.22)

    
    for i in range(1, h_data.GetNbinsX()+1):
        tmpVal = hCloneDivide.GetBinContent(i)
        if tmpVal != -0.:
            hDataDivide.SetBinContent(i, hDataDivide.GetBinContent(i)/tmpVal)
            hDataDivide.SetBinError(i, hDataDivide.GetBinError(i)/tmpVal)
            hCloneDivide.SetBinContent(i, hCloneDivide.GetBinContent(i)/tmpVal)
            hCloneDivide.SetBinError(i, hCloneDivide.GetBinError(i)/tmpVal)
            hDivide.SetBinContent(i, hDivide.GetBinContent(i)/tmpVal)
            hDivide.SetBinError(i, hDivide.GetBinError(i)/tmpVal)

            
    hCloneDivide.GetXaxis().SetTitleOffset(0.97)
    hCloneDivide.GetXaxis().SetLabelOffset(0.02)
    hCloneDivide.GetXaxis().SetTitle(xTitle)

    hCloneDivide.GetYaxis().SetNdivisions(504,rt.kTRUE)
    hCloneDivide.GetYaxis().SetTitleOffset(0.2)
    hCloneDivide.GetYaxis().SetTitleSize(0.22)
    hCloneDivide.GetYaxis().SetTitle("Data/Fit")
    hCloneDivide.GetXaxis().SetTicks("+")
    hCloneDivide.GetXaxis().SetTickLength(0.07)
    hCloneDivide.SetMarkerColor(rt.kBlue-10)
    hCloneDivide.Draw("e2")
    hDataDivide.Draw('pesame')
    hCloneDivide.Draw("axissame")


    pad2.Update()
    pad1.cd()
    pad1.Update()
    pad1.Draw()

    if tLeg==None:
        tLeg = rt.TLegend(0.7,0.6,0.9,0.8)
        tLeg.SetFillColor(0)
        tLeg.SetTextFont(42)
        tLeg.SetLineColor(0)
        tLeg.AddEntry(h_data,"Sim. Data","lep")
        tLeg.AddEntry(hClone,"Fit","lf")
    tLeg.Draw("same")

        
    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.05)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.15,0.9,"CMS simulation")
    l.DrawLatex(0.78,0.9,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.DrawLatex(0.2,0.8,boxLabel)

    c.cd()
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    
def print2DCanvas(c,h,printName):
    c.SetLogx(1)
    c.SetLogy(1)
    h.Draw("textcolz")
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    
def setFFColors(hNS):
    Red = array('d',  [0.00, 0.70, 0.90, 1.00, 1.00, 1.00, 1.00])
    Green = array('d',[0.00, 0.70, 0.90, 1.00, 0.90, 0.70, 0.00])
    Blue = array('d', [1.00, 1.00, 1.00, 1.00, 0.90, 0.70, 0.00])
    Length =array('d',[0.00, 0.20, 0.35, 0.50, 0.65, 0.8, 1.00]) # colors get darker faster at 4sigma
    rt.TColor.CreateGradientColorTable(7,Length,Red,Green,Blue,999)
    hNS.SetMaximum(5.1)
    hNS.SetMinimum(-5.1) # so the binning is 0 2 4
    hNS.SetContour(999)
    
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
    print stringMyStruct1+"};}"
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

def getBinSumDicts(x, y, z, btag=[[1,2,3]]):
    binXSumDict = {}
    binYSumDict = {}
    binXYSumDict = {}
    binXYZSumDict = {}
    
    for i in range(1,len(x)):
        binXSumDict[i] = ''
    for j in range(1,len(y)):
        binYSumDict[j] = ''
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            binXYSumDict[i,j] = ''
            
    iBinX = -1
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            for k in range(1,len(z)):
                iBinX += 1
                if (iBinX)%3+1 in btag: binXSumDict[i] += 'b%i+'%iBinX
                if (iBinX)%3+1 in btag: binYSumDict[j] += 'b%i+'%iBinX
                if (iBinX)%3+1 in btag: binXYSumDict[i,j] += 'b%i+'%iBinX
                binXYZSumDict[i,j,k] = 'b%i'%iBinX
                
    for i in range(1,len(x)):
        binXSumDict[i] = binXSumDict[i][0:-1]
    for j in range(1,len(y)):
        binYSumDict[j] = binYSumDict[j][0:-1]
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            binXYSumDict[i,j] = binXYSumDict[i,j][0:-1]
            
    return binXSumDict, binYSumDict, binXYSumDict, binXYZSumDict

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


def runToys(box,freq,w,outDir):
    
    setStyle()
    
    extRazorPdf = w.pdf('extRazorPdf')
    
    nll = extRazorPdf.createNLL(dataHist)
    m = rt.RooMinuit(nll)
    m.migrad()
    m.hesse()
    fr = m.save()
    #fr = extRazorPdf.fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Minimizer('Minuit2','improve'),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.SumW2Error(False))
    fr.Print('v')
    params = extRazorPdf.getParameters(dataHist)
    params.remove(w.var('BtagCut_TTj1b'))
    params.remove(w.var('BtagCut_TTj2b'))
    params.remove(w.var('BtagCut_TTj3b'))
    params.remove(w.var('MRCut_%s'%box))
    params.remove(w.var('RCut_%s'%box))
    
    hessePdf = fr.createHessePdf(params)
    minNll = rt.RooRealVar('minNll_%s'%box,'minNll_%s'%box,fr.minNll())
    minNll.setConstant(True)
    n2dll = rt.RooFormulaVar('n2dll_%s'%box,'2*@0-2*@1',rt.RooArgList(nll,minNll))
    
    rootTools.Utils.importToWS(w,n2dll)
    rootTools.Utils.importToWS(w,hessePdf)
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    
    th1x.setBins(nBins)


    for bkgd in ['TTj1b','TTj2b','TTj3b']:
        paramNames.remove('%s_%s_norm'%(box,bkgd))
        paramNames.append('Ntot_%s_%s'%(bkgd,box))
    paramNames.sort()

    unc = 'Bayes'
    if freq: unc = 'Freq'
    
    output = rt.TFile.Open(outDir+'/toys_%s_%s.root'%(unc,box),'recreate')
    output.cd()
    myTree = rt.TTree("myTree", "myTree")
    
    s1 = getTree(myTree,paramNames, nBins,box)
    
    value = setattr(s1, 'n2dll_%s'%box, n2dll.getVal())
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
    while iToy < nToys:
        if iToy%100==0: print "generating %ith toy"%iToy
        if freq:
            pSet = fr.floatParsFinal()
        elif bayes:
            pSet = fr.randomizePars()
        for p in rootTools.RootIterator.RootIterator(pSet):
            w.var(p.GetName()).setVal(p.getVal())
            w.var(p.GetName()).setError(p.getError())

        badPars = []
        for bkgd in ['TTj1b','TTj2b']:
            badPars.append(w.var('n_%s_%s'%(bkgd,box)).getVal() <= 0)
            badPars.append(w.var('b_%s_%s'%(bkgd,box)).getVal() <= 0)
            badPars.append(w.var('MR0_%s_%s'%(bkgd,box)).getVal() >= w.var('MR').getMin())
            badPars.append(w.var('R0_%s_%s'%(bkgd,box)).getVal()  >=  w.var('Rsq').getMin())
        if any(badPars): continue
        
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
        if bayes or freq:
            asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy_%i'%iToy))
            #asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('toy_%i'%iToy),rt.RooFit.Asimov())
        errorCountAfter = rt.RooMsgService.instance().errorCount()
        if errorCountAfter > errorCountBefore: continue

        pSetSave = pSet
        if freq:
            fr_toy = extRazorPdf.fitTo(asimov,rt.RooFit.Save(),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1),rt.RooFit.Warnings(False))
            if fr_toy.covQual() < 3: continue
            pSetSave = fr_toy.floatParsFinal()
            asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('pred_%i'%iToy),rt.RooFit.Asimov())
            #asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Name('pred_%i'%iToy))

        value = setattr(s1, 'n2dll_%s'%box, n2dll.getVal())

        for p in rootTools.RootIterator.RootIterator(pSetSave):
            value = setattr(s1, p.GetName(), p.getVal())

        for iBinX in range(0,nBins):
            th1x.setVal(iBinX+0.5)
            pdfValV = extRazorPdf.getValV(rt.RooArgSet(th1x)) * extRazorPdf.expectedEvents(rt.RooArgSet(th1x))
            if iToy%1000==0: print iBinX, asimov.weight(rt.RooArgSet(th1x)), pdfValV
            value = setattr(s1, 'b%i'%iBinX, float(asimov.weight(rt.RooArgSet(th1x))))
                        
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
    parser.add_option('-l','--lumi',dest="lumi", default=4000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-t','--toys',dest="nToys", default=10000,type="int",
                  help="number of toys")
    parser.add_option('--no-stat',dest="noStat",default=False,action='store_true',
                  help="just systematic uncertainty (no poisson fluctuations), default is statstical + systematic uncertainty")
    parser.add_option('--bayes',dest="bayes",default=True,action='store_true',
                  help="bayesian method")
    parser.add_option('--freq',dest="freq",default=False,action='store_true',
                  help="refit each toy with only statistical fluctuations, as in frequentist approach")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    
    box = options.box
    lumi = options.lumi
    nToys = options.nToys
    bayes = options.bayes
    freq = options.freq

    if freq: bayes = False
    else: bayes = True
    
    lumi_in = 0.
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get('w'+box)
            data = workspace.data('RMRTree')
            lumi_in = 1000.*float([g.replace('lumi-','') for g in f.split('_') if g.find('lumi')!=-1][0])

  
    w = rt.RooWorkspace("w"+box)
    
    paramNames, bkgs = initializeWorkspace(w,cfg,box)
    rootTools.Utils.importToWS(w,data)
    
    th1x = w.var('th1x')
    
    myTH1 = convertDataset2TH1(data, cfg, box, w)
    myTH1.Scale(lumi/lumi_in)
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), myTH1)
    rootTools.Utils.importToWS(w,dataHist)

    outputName = runToys(box,freq,w,options.outDir)
    output = rt.TFile.Open(outputName,"read")
    myTree = output.Get("myTree")

    c = rt.TCanvas('c','c',500,400)
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
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
        binXSumDict, binYSumDict, binXYSumDict, binXYZSumDict = getBinSumDicts(x,y,z,btag)
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
            
    setFFColors(h_FF_['1_2_3'])
    print2DCanvas(c,h_FF_['1_2_3'],options.outDir+'/h_FF.pdf')
    for btag in ['1','2','3']:
        setFFColors(h_FF_[btag])
        print2DCanvas(c,h_FF_[btag],'h_FF_'+btag)

    
    print1DCanvas(c,h_MR_['1_2_3'],h_data_MR_['1_2_3'],options.outDir+'/h_MR.pdf')
    print1DCanvas(c,h_Rsq_['1_2_3'],h_data_Rsq_['1_2_3'],options.outDir+'/h_Rsq.pdf')

    for btag in ['1','2','3']:
        print1DCanvas(c,h_MR_[btag],h_data_MR_[btag],options.outDir+'/h_MR_'+btag+'.pdf')
        print1DCanvas(c,h_Rsq_[btag],h_data_Rsq_[btag],options.outDir+'/h_Rsq_'+btag+'.pdf')

    
    unc = 'Bayes'
    if freq: unc = 'Freq'
        
    for h in [h_MR_['1_2_3'],h_Rsq_['1_2_3'],h_RsqMR_['1_2_3'],h_nBtagRsqMR]:
        h.SetDirectory(0)
    for h in [h_data_MR_['1_2_3'],h_data_Rsq_['1_2_3'],h_data_RsqMR_['1_2_3'],h_data_nBtagRsqMR]:
        h.SetDirectory(0)
        
    newoutput = rt.TFile(options.outDir+"/plots_"+unc+"_"+box+".root","recreate")
    newoutput.cd()
        
    for h in [h_MR_['1_2_3'],h_Rsq_['1_2_3'],h_RsqMR_['1_2_3'],h_nBtagRsqMR]:
        h.Write()
    for h in [h_data_MR_['1_2_3'],h_data_Rsq_['1_2_3'],h_data_RsqMR_['1_2_3'],h_data_nBtagRsqMR]:
        h.Write()
        
    output.Close()
    newoutput.Close()
    
    
