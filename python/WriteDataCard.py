from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import sys

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
    w.var('th1x').setBins(nBins)
    emptyHist3D = rt.TH3D("emptyHist3D","emptyHist3D",len(x)-1,x,len(y)-1,y,len(z)-1,z)

    commands = cfg.getVariables(box, "combine_pdfs")
    for command in commands:
        if command.find('SUM::')!=-1:
            w.factory(command)
        else:
            myclass = command.split('::')[0]
            remaining = command.split('::')[1]
            name = remaining.split('(')[0]
            mytuple = remaining.replace(name,'').replace('(','').replace(')','')
            mylist = mytuple.split(',')
            arglist = [name, name]
            for myvar in mylist:
                arglist.append(w.var(myvar))
            args = tuple(arglist)
            pdf = getattr(rt,myclass)(*args)
            if hasattr(pdf,'setTH3Binning'):
                pdf.setTH3Binning(emptyHist3D)
            rootTools.Utils.importToWS(w,pdf)

    return paramNames


def initializeWorkspace_noFit(w,cfg,box):
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
    data.fillHistogram(myTH3, varList,"MR>%f && MR<%f && Rsq>%f && Rsq<%f && nBtag >= %f && nBtag <= %f"%(x[0],x[-1],y[0],y[-1],z[0],z[-1]))
    data.fillHistogram(myTH2, varList2D,"MR>%f && MR<%f && Rsq>%f && Rsq<%f && nBtag >= %f && nBtag <= %f"%(x[0],x[-1],y[0],y[-1],z[0],z[-1]))
    
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

        
def writeDataCard_noFit(box,model,txtfileName,bkgs,paramNames,w):
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
                        (w.data("%s_%s"%(box,model)).sumEntries(),w.data("%s_%s"%(box,"TTj1b")).sumEntries(),
                        w.data("%s_%s"%(box,"TTj2b")).sumEntries(), w.data("%s_%s"%(box,"TTj3b")).sumEntries()))
        
        txtfile.write("------------------------------------------------------------\n")
        txtfile.write("lumi			lnN	%.3f       1.00	1.00 1.00\n"%(1.05))
        txtfile.write("ttj1b_%s			lnN	1.00       %.3f	1.00 1.00\n"%(box,1.1))
        txtfile.write("ttj2b_%s			lnN	1.00       1.00	%.3f 1.00\n"%(box,1.1))
        txtfile.write("ttj3b_%s			lnN	1.00       1.00	1.00 %.3f\n"%(box,1.1))
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
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (use MC directly)")
    parser.add_option('--asimov',dest="asimov",default=False,action='store_true',
                  help="use asimov dataset derived from fit")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box = options.box
    lumi = options.lumi
    noFit = options.noFit

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
    
    if noFit:
        paramNames = initializeWorkspace_noFit(w,cfg,box)
    else:
        paramNames = initializeWorkspace(w,cfg,box)
    
    rootTools.Utils.importToWS(w,data)
    
    th1x = w.var('th1x')
    

    myTH1 = convertDataset2TH1(data, cfg, box, w)
    myTH1.Scale(lumi/lumi_in)
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), myTH1)
    if not options.asimov:
        rootTools.Utils.importToWS(w,dataHist)
    
    if noFit:
        data1b = data.reduce("nBtag==1")
        myTH11b = convertDataset2TH1(data1b, cfg, box, w)
        myTH11b.Scale(lumi/lumi_in)
        dataHist1b = rt.RooDataHist("%s_%s"%(box,"TTj1b"),"%s_%s"%(box,"TTj1b"),rt.RooArgList(th1x), myTH11b)
        rootTools.Utils.importToWS(w,dataHist1b)
        
        data2b = data.reduce("nBtag==2")
        myTH12b = convertDataset2TH1(data2b, cfg, box, w)
        myTH12b.Scale(lumi/lumi_in)
        dataHist2b = rt.RooDataHist("%s_%s"%(box,"TTj2b"),"%s_%s"%(box,"TTj2b"),rt.RooArgList(th1x), myTH12b)
        rootTools.Utils.importToWS(w,dataHist2b)
        
        data3b = data.reduce("nBtag==3")
        myTH13b = convertDataset2TH1(data3b, cfg, box, w)
        myTH13b.Scale(lumi/lumi_in)
        dataHist3b = rt.RooDataHist("%s_%s"%(box,"TTj3b"),"%s_%s"%(box,"TTj3b"),rt.RooArgList(th1x), myTH13b)
        rootTools.Utils.importToWS(w,dataHist3b)

    elif options.asimov:
        ntot_ttj1b = rt.RooRealVar('Ntot_TTj1b','Ntot_TTj1b',data.sumEntries("nBtag==1"),0,10000)
        ntot_ttj2b = rt.RooRealVar('Ntot_TTj2b','Ntot_TTj2b',data.sumEntries("nBtag==2"),0,10000)
        ntot_ttj3b = rt.RooRealVar('Ntot_TTj3b','Ntot_TTj3b',data.sumEntries("nBtag==3"),0,10000)
        coefList = rt.RooArgList()
        coefList.add(ntot_ttj1b)
        coefList.add(ntot_ttj2b)
        coefList.add(ntot_ttj3b)
        pdfList = rt.RooArgList()
        pdfList.add(w.pdf('%s_%s'%(box,'TTj1b')))
        pdfList.add(w.pdf('%s_%s'%(box,'TTj2b')))
        pdfList.add(w.pdf('%s_%s'%(box,'TTj3b')))
        extRazorPdf = rt.RooAddPdf('extRazorPdf','extRazorPdf',pdfList,coefList)
        fr = extRazorPdf.fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Minimizer('Miniuit2','migrad'),rt.RooFit.PrintLevel(-1),rt.RooFit.PrintEvalErrors(-1))
        fr.Print('v')
        asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Asimov(),rt.RooFit.Name('data_obs'))
        rootTools.Utils.importToWS(w,asimov)
    
    
    sigTH1 = convertDataset2TH1(signalDs, cfg, box, w,"signal")
    sigTH1.Scale(lumi/lumi_in)
    sigDataHist = rt.RooDataHist('%s_%s'%(box,model),'%s_%s'%(box,model),rt.RooArgList(th1x), sigTH1)
    rootTools.Utils.importToWS(w,sigDataHist)

    w.Print('v')
            
    outFile = 'razor_combine_%s_%s_lumi-%.1f_%s.root'%(model,massPoint,lumi/1000.,box)
    
    outputFile = rt.TFile.Open(options.outDir+"/"+outFile,"recreate")
    w.Write()
    if noFit:
        writeDataCard_noFit(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),["TTj1b","TTj2b","TTj3b"],paramNames,w)
    else:
        writeDataCard(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),["TTj1b","TTj2b","TTj3b"],paramNames,w)
    os.system("cat %s"%options.outDir+"/"+outFile.replace(".root",".txt"))
