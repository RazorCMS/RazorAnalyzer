from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import sys
from itertools import *
from operator import *

seed = 1988


def fixPars(w, label, doFix=True, setVal=None):
    parSet = w.allVars()
    for par in rootTools.RootIterator.RootIterator(parSet):
        if label in par.GetName():
            par.setConstant(doFix)
            if setVal is not None: par.setVal(setVal)

def initializeWorkspace(w,cfg,box,scaleFactor=1.,x=None,y=None,z=None):
    parameters = cfg.getVariables(box, "combine_parameters")
    paramNames = []
    for parameter in parameters:
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if not ("Cut" in paramName or "Ntot" in paramName):
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
            
        # fix Rsq MR cut parameters
        fixPars(w,"Cut")

        # float normalization parameters
        fixPars(w,"Ntot",False)
        
        # turn off shape parameters if no events in b-tag bin:
        for i in [0, 1, 2, 3]:
            if "Ntot" in paramName and "%ib"%i in paramName:
                w.var(paramName).setVal(scaleFactor * (w.data("RMRTree").sumEntries("nBtag==%i"%i) ))
                if not w.var(paramName).getVal():
                    fixPars(w,"%ib"%i)    
    if x is None or y is None or z is None:
        x = array('d', cfg.getBinning(box)[0]) # MR binning
        y = array('d', cfg.getBinning(box)[1]) # Rsq binning
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)

    w.factory('th1x[0,0,%i]'%nBins)
    w.var('th1x').setBins(nBins)
    emptyHist3D = rt.TH3D("emptyHist3D","emptyHist3D",len(x)-1,x,len(y)-1,y,len(z)-1,z)

    iBinX = -1
    sidebandBins = []
    for ix in range(1,len(x)):
        for iy in range(1,len(y)):
            for iz in range(1,len(z)):
                iBinX+=1
                emptyHist3D.SetBinContent(ix,iy,iz,1.)
                w.var('MR').setVal(emptyHist3D.GetXaxis().GetBinCenter(ix))
                w.var('Rsq').setVal(emptyHist3D.GetYaxis().GetBinCenter(iy))
                w.var('nBtag').setVal(emptyHist3D.GetZaxis().GetBinCenter(iz))
                inSideband = ( w.var('MR').inRange('LowMR') * w.var('Rsq').inRange('LowMR') * w.var('nBtag').inRange('LowMR') )
                inSideband += ( w.var('MR').inRange('LowRsq') * w.var('Rsq').inRange('LowRsq') * w.var('nBtag').inRange('LowRsq') )
                if inSideband: sidebandBins.append(iBinX)
                    
    sidebandGroups = []
    for k, g in groupby(enumerate(sidebandBins), lambda (i,x):i-x):
        consecutiveBins = map(itemgetter(1), g)
        sidebandGroups.append([consecutiveBins[0],consecutiveBins[-1]+1])

    #for iSideband, sidebandGroup in enumerate(sidebandGroups):
    #    w.var("th1x").setRange("sideband%i"%iSideband,sidebandGroup[0],sidebandGroup[1])        

    w.Print('v')
    commands = cfg.getVariables(box, "combine_pdfs")
    bkgs = []
    for command in commands:
        lower = command.lower()
        if lower.find('sum::')!=-1 or lower.find('prod::')!=-1 or lower.find('expr::')!=-1:
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
            bkg = name.split("_")
            if box in bkg: bkg.remove(box)
            bkgs.append("_".join(bkg))
    return paramNames, bkgs


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

def convertDataset2TH1(data, cfg, box, workspace, th1Name = 'h', x = array('d',[]), y = array('d',[]), z = array('d',[])):
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

    if cfg is not None:    
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
        obsRate = w.data("data_obs").sumEntries()
        nBkgd = len(bkgs)
        rootFileName = txtfileName.replace('.txt','.root')
        rates = [w.data("%s_%s"%(box,model)).sumEntries()]
        rates.extend([w.var('Ntot_%s_%s'%(bkg,box)).getVal() for bkg in bkgs])
        processes = ["%s_%s"%(box,model)]
        processes.extend(["%s_%s"%(box,bkg) for bkg in bkgs])
        lumiErrs = [1.05]
        lumiErrs.extend([1 for bkg in bkgs])
        divider = "------------------------------------------------------------\n"
        datacard = "imax 1 number of channels\n" + \
                   "jmax %i number of backgrounds\n"%nBkgd + \
                   "kmax * number of nuisance parameters\n" + \
                   divider + \
                   "observation	%.3f\n"%obsRate + \
                   divider + \
                   "shapes * * %s w%s:$PROCESS w%s:$PROCESS_$SYSTEMATIC\n"%(rootFileName,box,box) + \
                   divider
        binString = "bin"
        processString = "process"
        processNumberString = "process"
        rateString = "rate"
        lumiString = "lumi\tlnN"
        for i in range(0,len(bkgs)+1):
            binString +="\t%s"%box
            processString += "\t%s"%processes[i]
            processNumberString += "\t%i"%i
            rateString += "\t%.3f" %rates[i]
            lumiString += "\t%.3f"%lumiErrs[i]
        binString+="\n"; processString+="\n"; processNumberString+="\n"; rateString +="\n"; lumiString+="\n"
        datacard+=binString+processString+processNumberString+rateString+divider
        # now nuisances
        datacard+=lumiString
        for paramName in paramNames:
            datacard += "%s  	flatParam\n"%(paramName)
        txtfile = open(txtfileName,"w")
        txtfile.write(datacard)
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
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('--no-fit',dest="noFit",default=False,action='store_true',
                  help="Turn off fit (use MC directly)")
    parser.add_option('--fit',dest="fit",default=False,action='store_true',
                  help="perform a fit first")
    parser.add_option('--print-yields',dest="printYields",default=False,action='store_true',
                  help="print yields")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box = options.box
    lumi = options.lumi
    noFit = options.noFit
    printYields = options.printYields

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
    
    rootTools.Utils.importToWS(w,data)
    
    if noFit:
        paramNames = initializeWorkspace_noFit(w,cfg,box)
    else:
        paramNames, bkgs = initializeWorkspace(w,cfg,box,lumi/lumi_in)
    
    
    th1x = w.var('th1x')
    
    myTH1 = convertDataset2TH1(data, cfg, box, w)
    myTH1.Scale(lumi/lumi_in)
    dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), rt.RooFit.Import(myTH1))
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

    elif options.fit:
        fr = w.pdf('extRazorPdf').fitTo(dataHist,rt.RooFit.Save(),rt.RooFit.Minimizer('Minuit2','migrad'),rt.RooFit.PrintLevel(-1),rt.RooFit.SumW2Error(False),rt.RooFit.PrintEvalErrors(-1))
        fr.Print('v')
        #asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),rt.RooFit.Asimov(),rt.RooFit.Name('data_obs'))
        #rootTools.Utils.importToWS(w,asimov)
    
    
    sigTH1 = convertDataset2TH1(signalDs, cfg, box, w,"signal")
    sigTH1.Scale(lumi/lumi_in)
    sigDataHist = rt.RooDataHist('%s_%s'%(box,model),'%s_%s'%(box,model),rt.RooArgList(th1x), sigTH1)
    rootTools.Utils.importToWS(w,sigDataHist)
    
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    btagMin = z[0]
    btagMax = z[-1]        
    if btagMax>btagMin+1:            
        outFile = 'razor_combine_%s_%s_lumi-%.1f_%i-%ibtag_%s.root'%(model,massPoint,lumi/1000.,btagMin,btagMax-1,box)
    else:
        outFile = 'razor_combine_%s_%s_lumi-%.1f_%ibtag_%s.root'%(model,massPoint,lumi/1000.,btagMin,box)
    
    outputFile = rt.TFile.Open(options.outDir+"/"+outFile,"recreate")
    w.Write()
    if noFit:
        writeDataCard_noFit(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),["TTj1b","TTj2b","TTj3b"],paramNames,w)
    else:
        writeDataCard(box,model,options.outDir+"/"+outFile.replace(".root",".txt"),bkgs,paramNames,w)
    os.system("cat %s"%options.outDir+"/"+outFile.replace(".root",".txt"))

    

    if printYields:
        x = array('d', cfg.getBinning(box)[0]) # MR binning
        y = array('d', cfg.getBinning(box)[1]) # Rsq binning
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
        zBins = len(z)-1
        yBins = len(y)-1
        xBins = len(x)-1
        csvFile = open(options.outDir+"/"+outFile.replace(".root",".csv"),'w')
        csvOutput = "bin number,MR range,Rsq range,b-tags,signal yield (S),background yield from MC (B),background yield from fit (F),S/B,S/F"
        for iBin in range(0,th1x.getBins()):
            signal = sigDataHist.sumEntries("th1x>=%i && th1x<%i+1"%(iBin,iBin))
            mc = w.data("data_obs").sumEntries("th1x>=%i && th1x<%i+1"%(iBin,iBin))
            asimov = w.pdf("extRazorPdf").generateBinned(rt.RooArgSet(th1x),rt.RooFit.Asimov())
            fit = asimov.sumEntries("th1x>=%i && th1x<%i+1"%(iBin,iBin))
            zBin = iBin % zBins
            yBin = ( (iBin - zBin)/(zBins) ) % (yBins)
            xBin =  (iBin - zBin - yBin*zBins ) / (zBins*yBins)
            csvOutput += "\n%i,%i-%i,%.2f-%.2f,%i,%f,%f,%f,=E%i/F%i,=E%i/G%i" % (iBin, x[xBin],x[xBin+1],y[yBin],y[yBin+1],z[zBin],signal,mc,fit,iBin+2,iBin+2,iBin+2,iBin+2)
        csvFile.write(csvOutput)
        csvFile.close()
        print ""
        print "yields written into %s"%(options.outDir+"/"+outFile.replace(".root",".csv"))
