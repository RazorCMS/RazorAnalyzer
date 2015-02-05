from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os

seed = 1988

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

def convertDataset2TH1(data, cfg, box, workspace, th1Name = 'h'):
    """Get the cocktail dataset from the file"""
    
    row = data.get()

    MR = row['MR']
    Rsq = row['Rsq']
    nBtag = row['nBtag']
    
    varSet = rt.RooArgSet(MR,Rsq,nBtag)
    varList = rt.RooArgList(MR,Rsq,nBtag)
    varList2D = rt.RooArgList(MR,Rsq)
        
    mRmin = row['MR'].getMin()
    mRmax = row['MR'].getMax()
    rsqMin = row['Rsq'].getMin()
    rsqMax = row['Rsq'].getMax()
    nbtagMin = row['nBtag'].getMin()
    nbtagMax = row['nBtag'].getMax()

    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    myTH3 = rt.TH3D(th1Name+box, th1Name+box, len(x)-1, x, len(y)-1, y, len(z)-1, z)
    myTH2 = rt.TH2D(th1Name+box+"2d", th1Name+box+"2d", len(x)-1, x, len(y)-1, y)
    myTH2.Sumw2()

    # fills automatically with weight
    data.fillHistogram(myTH3, varList,"MR>0")
    data.fillHistogram(myTH2, varList2D,"MR>0")
    
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    myTH1 = rt.TH1D(th1Name+box+"1d",th1Name+box+"1d",nBins,0,nBins)
    i = 0
    for ix in range(1,len(x)):
        for iy in range(1,len(y)):
            for iz in range(1,len(z)):
                i+= 1
                myTH1.SetBinContent(i,myTH3.GetBinContent(ix,iy,iz))
   
    return myTH1


def writeDataCard(box,model,txtfileName,bkgs,paramNames,w):
        txtfile = open(txtfileName,"w")
        txtfile.write("imax 1 number of channels\n")
        nBkgd = 3
        txtfile.write("jmax %i number of backgrounds\n"%nBkgd)
        txtfile.write("kmax * number of nuisance parameters\n")
        txtfile.write("------------------------------------------------------------\n")
        txtfile.write("observation	%.3f\n"%
                      w.data("data_obs").sumEntries())
        txtfile.write("------------------------------------------------------------\n")
        txtfile.write("shapes * * %s w%s:$PROCESS w%s:$PROCESS_$SYSTEMATIC\n"%
                      (txtfileName.replace('.txt','.root'),box,box))
        txtfile.write("------------------------------------------------------------\n")
        txtfile.write("bin		%s			%s			%s			%s\n"%(box,box,box,box))
        txtfile.write("process		%s_%s 	%s_%s	%s_%s	%s_%s\n"%
                        (box,model,box,bkgs[0],box,bkgs[1],box,bkgs[2]))
        txtfile.write("process        	0          		1			2			3\n")
        txtfile.write("rate            %.3f		%.3f		%.3f		%.3f\n"%
                        (w.data("%s_%s"%(box,model)).sumEntries(),w.data("RMRTree").sumEntries("nBtag==1")*lumi/lumi_in,
                        w.data("RMRTree").sumEntries("nBtag==2")*lumi/lumi_in,w.data("RMRTree").sumEntries("nBtag==3")*lumi/lumi_in))
        
        txtfile.write("------------------------------------------------------------\n")
        txtfile.write("lumi			lnN	%.3f       1.00	1.00 1.00\n"%(1.05))
        for paramName in paramNames:
                txtfile.write("%s  	flatParam\n"%
                              (paramName))
        txtfile.close()
        
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

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box = options.box
    lumi = options.lumi

    lumi_in = 0.
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get('w'+box)
            if f.lower().find('t1')!=-1 or f.lower().find('t2')!=-1:
                signalDs = workspace.data('RMRTree')
                model = f.split('-')[1].split('_')[0]
                massPoint = '_'.join(f.split('_')[3:5])
            else:
                data = workspace.data('RMRTree')
            lumi_in = 1000.*float([g.replace('lumi-','') for g in f.split('_') if g.find('lumi')!=-1][0])

  
    w = rt.RooWorkspace("w"+box)
    paramNames = initializeWorkspace(w,cfg,box)
    
    rootTools.Utils.importToWS(w,data)
    
    th1x = w.var('th1x')
    
    myTH1 = convertDataset2TH1(data, cfg, box, w)
    myTH1.Scale(lumi/lumi_in)
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), myTH1)
    rootTools.Utils.importToWS(w,dataHist)

    sigTH1 = convertDataset2TH1(signalDs, cfg, box, w,"signal")
    sigTH1.Scale(lumi/lumi_in)
    sigDataHist = rt.RooDataHist('%s_%s'%(box,model),'%s_%s'%(box,model),rt.RooArgList(th1x), sigTH1)
    rootTools.Utils.importToWS(w,sigDataHist)


    w.Print('v')

            
    outFile = 'razor_combine_%s_%s_lumi-%.1f_%s.root'%(model,massPoint,lumi/1000.,box)
        
    
    outputFile = rt.TFile.Open(options.outDir+"/"+outFile,"recreate")
    w.Write()
    writeDataCard(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),["TTj1b","TTj2b","TTj3b"],paramNames,w)
    os.system("cat %s"%options.outDir+"/"+outFile.replace(".root",".txt"))