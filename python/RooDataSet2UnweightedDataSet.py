from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *

seed = 1988

def initializeWorkspace(w,cfg):
    variables = cfg.getVariablesRange(box,"variables",w)
    parameters = cfg.getVariables(box, "parameters")
    paramNames = []
    for parameter in parameters:
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if paramName.find("Cut")==-1 and paramName.find("Ntot")==-1:
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
        else:
            if paramName.find("Ntot")==-1:
                w.var(paramName).setConstant(True)
            else:
                w.var(paramName).setConstant(False)
                
def convertDataset2UnweightedToy(data, cfg, box, workspace, uwName = 'uw'):
    """Get the cocktail dataset from the file"""
    row = data.get()

    MR = row['MR']
    Rsq = row['Rsq']
    nBtag = row['nBtag']
    
    varSet = rt.RooArgSet(MR,Rsq,nBtag)
    varList = rt.RooArgList(MR,Rsq,nBtag)
    varList2D = rt.RooArgList(MR,Rsq)
    uwdata = rt.RooDataSet(uwName+'tree','Unweighted Cocktail',varSet)
        
    mRmin = row['MR'].getMin()
    mRmax = row['MR'].getMax()
    rsqMin = row['Rsq'].getMin()
    rsqMax = row['Rsq'].getMax()
    nbtagMin = row['nBtag'].getMin()
    nbtagMax = row['nBtag'].getMax()
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    myTH3 = rt.TH3D(uwName+box, uwName+box, 100, mRmin, mRmax, 70, rsqMin, rsqMax, 3, nbtagMin, nbtagMax)
    myTH2 = rt.TH2D(uwName+box+"2d", uwName+box+"2d", 100, mRmin, mRmax, 70, rsqMin, rsqMax)
    myTH2.Sumw2()

    # fills automatically with weight
    data.fillHistogram(myTH3, varList,"MR>0")
    data.fillHistogram(myTH2, varList2D,"MR>0")
    
    c = rt.TCanvas("c","c",600,400)
    #rt.gStyle.SetOptStat(1001000011)
    rt.gStyle.SetOptStat(0)
    myTH2.SetTitle("Weighted %s"%box)
    sumW2 = 0
    for i in range(0,wdata.numEntries()):
       wdata.get(i)
       sumW2+=(wdata.weight())*(wdata.weight())

    print "sum (weights)^2 = %.1f" %sumW2
    print "(sum weights)^2 = %.1f" %((wdata.sumEntries())*(wdata.sumEntries()))
    effEntries = (((wdata.sumEntries())*(wdata.sumEntries()))/sumW2)
    print "effective entries = %.1f"%effEntries
    myTH2.GetXaxis().SetTitle("M_{R}")
    myTH2.GetYaxis().SetTitle("R^{2}")
    myTH2.GetXaxis().SetMoreLogLabels()
    myTH2.GetYaxis().SetMoreLogLabels()
    myTH2.GetXaxis().SetNoExponent()
    myTH2.GetYaxis().SetNoExponent()
    myTH2.Draw("colz")
    c.SetLogy()
    c.SetLogx()

    
    #l = rt.TLatex()
    #l.SetTextAlign(11)
    #l.SetTextSize(0.045)
    #l.SetTextFont(42)
    #l.SetNDC()
    #l.DrawLatex(0.18,0.84,"CMS Simulation 5 fb  ^{-1} (13 TeV)")
    #l.DrawLatex(0.18,0.77,"Razor MultiJet Box, QCD")
    #l.DrawLatex(0.11,0.84,"CMS Simulation 5 fb  ^{-1} (13 TeV)")
    #l.DrawLatex(0.11,0.77,"Razor MultiJet Box")
    #l.DrawLatex(0.56,0.84,"pp #rightarrow #tilde{g}#tilde{g},  #tilde{g}#rightarrowb#bar{b}#tilde{#chi}^{0}_{1}")
    #l.DrawLatex(0.5,0.77,"m_{#tilde{g}} = %i GeV, m  _{#tilde{#chi}} = %i GeV"%(1000,900))
    
    c.Print(options.outDir+"/TH2D_SMCocktail_weighted_%s.pdf"%box)
    c.Print(options.outDir+"/TH2D_SMCocktail_weighted_%s.C"%box)
    
    print wdata.weight()
    Nev = myTH3.Integral()
    Nent = myTH3.GetEntries()
    print "weighted events %.1f"% Nev
    print "entries  %d"% Nent
    Npois = rt.RooRandom.randomGenerator().Poisson(Nev)
    for i in range(0,Npois):
       myMR = rt.Double()
       myRsq = rt.Double()
       mynBtag = rt.Double()
       myTH3.GetRandom3(myMR,myRsq,mynBtag)
       mynBtag = int(mynBtag)
       varSet.setRealValue('MR',myMR)
       varSet.setRealValue('Rsq',myRsq)
       varSet.setRealValue('nBtag',mynBtag)
       uwdata.add(varSet)
    

    myTH2Toy = rt.TH2D("h", "h", 100, mRmin, mRmax, 70, rsqMin, rsqMax)
    uwdata.fillHistogram(myTH2Toy, varList2D,"MR>0")
    myTH2Toy.SetTitle("Unweighted %s"%box)
    myTH2Toy.GetXaxis().SetTitle("M_{R}")
    myTH2Toy.GetYaxis().SetTitle("R^{2}")
    myTH2Toy.GetXaxis().SetMoreLogLabels()
    myTH2Toy.GetYaxis().SetMoreLogLabels()
    myTH2Toy.GetXaxis().SetNoExponent()
    myTH2Toy.GetYaxis().SetNoExponent()
    myTH2Toy.Draw("colz")
    c.Print(options.outDir+"/TH2D_SMCocktail_unweighted_%s.pdf"%box)
    c.Print(options.outDir+"/TH2D_SMCocktail_unweighted_%s.C"%box)

    return uwdata


    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box

    ids = 0
    
    w = rt.RooWorkspace("w"+box)

    initializeWorkspace(w,cfg)
    
    ds = []
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get('w'+box)
            ds.append(workspace.data('RMRTree').Clone('RMRTree_%i'%ids))
            ids+=1

    
    wdata = ds[0].Clone('RMRTree')
    for ids in range(1,len(ds)):
        wdata.append(ds[ids])
    
    uwdata = convertDataset2UnweightedToy(wdata, cfg, box, w, uwName = 'uw')

    uwdata.SetName('RMRTree')
    rootTools.Utils.importToWS(w,uwdata)

    
    inFiles = [f for f in args if f.lower().endswith('.root')]
            
    if len(inFiles)==1:
        outFile = inFiles[0].split('/')[-1].replace('weighted','unweighted')
        
    outFile = rt.TFile.Open(options.outDir+"/"+outFile,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
    
            
