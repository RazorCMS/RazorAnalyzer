from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
import sys
from array import *

k_T = 689.1/424.5
k_Z = 3.*2008.4/5482.
k_W = 3.*20508.9/50100.0

k_QCD = {}
    
boxes = {'MuEle':'(box==0)',
         'MuMu':'(box==1)',
         'EleEle':'(box==2)',
         'MuSixJet':'(box==3)',
         'MuFourJet':'(box==4)',
         'MuMultiJet':'(box==3||box==4)',
         'MuJet':'(box==5)',
         'EleSixJet':'(box==6)',
         'EleFourJet':'(box==7)',
         'EleMultiJet':'(box==6||box==7)',
         'EleJet':'(box==8)',
         'LooseLeptonSixJet':'(box==9)',
         'LooseLeptonFourJet':'(box==10)',
         'LooseLeptonMultiJet':'(box==9||box==10)',
         'SixJet':'(box==11)',
         'FourJet':'(box==12)',
         'FourToSixJet':'((box==11||box==12)&&(nSelectedJets>=4&&nSelectedJets<7))',
         'SevenJet':'((box==11||box==12)&&(nSelectedJets>=7))',
         'MultiJet':'(box==11||box==12)',
         'LooseLeptonDiJet':'(box==13)',
         'DiJet':'(box==14)'}

backgrounds = ['dyjetstoll_htbinned', 'qcd_htbinned', 'ttjets', 'zjetstonunu_htbinned', 'multiboson', 'singletop', 'wjetstolnu_htbinned']

dPhiCut = 2.7

MTCut = -1

def initializeWorkspace(w,cfg,box):
    variables = cfg.getVariablesRange(box,"variables",w)
    
    w.factory('W[1.,0.,+INF]')
    w.set('variables').add(w.var('W'))
    return w

def getSumOfWeights(tree, cfg, box, workspace, useWeight, f, lumi, lumi_in):
    if f.find('SMS')!=-1:
        k = 1.
    elif f.find('TTJets')!=-1:
        k = k_T
    elif f.find('DYJets')!=-1 or f.find('ZJets')!=-1:
        k = k_Z
    elif f.find('WJets')!=-1:
        k = k_W
    else:
        k = 1.
        
    args = workspace.set("variables")
    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin()
    mRmax = args['MR'].getMax()

    #we cut away events outside our Rsq window
    rsqMin = args['Rsq'].getMin()
    rsqMax = args['Rsq'].getMax()

    btagMin =  args['nBtag'].getMin()
    btagMax =  args['nBtag'].getMax()
    
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    btagCutoff = 3
    if box in ["MuEle", "MuMu", "EleEle"]:
        btagCutoff = 1
        
    boxCut = boxes[box]

    label = f.replace('.root','').split('/')[-1]
    htemp = rt.TH1D('htemp_%s'%label,'htemp_%s'%label,len(z)-1,z)

    if useWeight:
        tree.Project(htemp.GetName(),
                    'min(nBTaggedJets,%i)'%btagCutoff,
                    '(%f/%f) * %f * weight * (MR > %f && MR < %f && Rsq > %f && Rsq < %f && min(nBTaggedJets,%i) >= %i && min(nBTaggedJets,%i) < %f && %s && abs(dPhiRazor) < %f)' % (lumi,lumi_in,k,mRmin,mRmax,rsqMin,rsqMax,btagCutoff,btagMin,btagCutoff,btagMax,boxCut,dPhiCut))
    else:
        tree.Project(htemp.GetName(),
                    'MR',
                    '(MR > %f && MR < %f && Rsq > %f && Rsq < %f && min(nBTaggedJets,%i) >= %i && min(nBTaggedJets,%i) < %f && %s && abs(dPhiRazor) < %f)' % (mRmin,mRmax,rsqMin,rsqMax,btagCutoff,btagMin,btagCutoff,btagMax,boxCut,dPhiCut))
        
    return [htemp.GetBinContent(i) for i in range(1,len(z))]
        
    
def convertTree2TH1(tree, cfg, box, workspace, useWeight, f, lumi, lumi_in, treeName):
    """Create 3D histogram for direct use with Combine"""
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    #get k factor for each btag bin, adjusted to correct for QCD
    if 'SMS' in f:
        k = [1. for z_bin in z[:-1]]
    elif 'TTJets' in f:
        k = [k_T*k_btag for k_btag in k_QCD[box]]
    elif 'DYJets' in f or 'ZJets' in f:
        k = [k_Z*k_btag for k_btag in k_QCD[box]]
    elif 'WJets' in f:
        k = [k_W*k_btag for k_btag in k_QCD[box]]
    else:
        k = k_QCD[box]

    #get variables and limits
    args = workspace.set("variables")
    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin()
    mRmax = args['MR'].getMax()

    #we cut away events outside our Rsq window
    rsqMin = args['Rsq'].getMin()
    rsqMax = args['Rsq'].getMax()

    btagMin =  args['nBtag'].getMin()
    btagMax =  args['nBtag'].getMax()
    
    #make histogram with razor binning
    label = f.replace('.root','').split('/')[-1]
    myTH3 = rt.TH3D(treeName+"3d",treeName+"3d",len(x)-1,x,len(y)-1,y,len(z)-1,z)

    #temp histogram for btag bins
    htemp = rt.TH1D('htemp_%s'%label,'htemp_%s'%label,len(z)-1,z)

    btagCutoff = 3
    if box in ["MuEle", "MuMu", "EleEle"]:
        btagCutoff = 1
        
    boxCut = boxes[box]
    
    #get list of entries passing the cuts
    tree.Draw('>>elist',
              'MR > %f && MR < %f && Rsq > %f && Rsq < %f && min(nBTaggedJets,%i) >= %i && min(nBTaggedJets,%i) < %i && %s && abs(dPhiRazor) < %f' % (mRmin,mRmax,rsqMin,rsqMax,btagCutoff,btagMin,btagCutoff,btagMax,boxCut,dPhiCut),
              'entrylist')
        
    elist = rt.gDirectory.Get('elist')
    
    #loop and fill histogram
    entry = -1;
    numEntriesByBtag = [0 for i in range(len(z)-1)]
    sumEntriesByBtag = [0. for i in range(len(z)-1)]
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)

        #get weight and fill
        nBTags = min(tree.nBTaggedJets,btagCutoff)
        btag_bin = htemp.FindBin(nBTags) - 1
        if useWeight:
            theWeight = tree.weight*lumi*k[btag_bin]/lumi_in
            myTH3.Fill(tree.MR, tree.Rsq, nBTags, theWeight)
            numEntriesByBtag[btag_bin] += 1
            sumEntriesByBtag[btag_bin] += theWeight
        else:
            myTH3.Fill(tree.MR, tree.Rsq, nBTags)
            numEntriesByBtag[btag_bin] += 1

    #unroll into TH1D
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    myTH1 = rt.TH1D(treeName,treeName,nBins,0,nBins)
    i = 0
    for ix in range(1,len(x)):
        for iy in range(1,len(y)):
            for iz in range(1,len(z)):
                i+= 1
                myTH1.SetBinContent(i,myTH3.GetBinContent(ix,iy,iz))

    print "Filename: %s"%f
    print "Sample: %s"%treeName
    print "Scale Factors     [ %s ] ="%box, k
    print "Number of Entries [ %s ] ="%(box), numEntriesByBtag
    print "Sum of Weights    [ %s ] ="%(box), sumEntriesByBtag

    myTH1.SetDirectory(0) #prevent it from going out of scope
    return myTH1

def writeDataCard_th1(box,model,txtfileName,hists):
    bkgs = [bkg for bkg in backgrounds if bkg in hists]
    obsRate = hists["data_obs"].Integral()
    nBkgd = len(bkgs)
    rootFileName = txtfileName.replace('.txt','.root')
    rates = [hists[model].Integral()]
    rates.extend([hists[bkg].Integral() for bkg in bkgs])
    processes = [model]
    processes.extend(bkgs)
    lumiErrs = [1.05]
    lumiErrs.extend([1.05 for bkg in bkgs])
    mcErrs = {} #dictionary of uncorrelated mc bkgd lnN uncertainties
    for bkg in bkgs:
            mcErrs[bkg] = [1.00]
            mcErrs[bkg].extend([1.00 + 0.10*(bkg==bkg1) for bkg1 in bkgs])
            
    divider = "------------------------------------------------------------\n"
    datacard = "imax 1 number of channels\n" + \
               "jmax %i number of backgrounds\n"%nBkgd + \
               "kmax * number of nuisance parameters\n" + \
               divider + \
               "observation	%.3f\n"%obsRate + \
               divider + \
               "shapes * * %s $PROCESS $PROCESS_$SYSTEMATIC\n"%(rootFileName) + \
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
        
    mcErrStrings = {}
    for bkg in bkgs:
            mcErrStrings[bkg] = "%s_norm\tlnN"%(bkg)
            for i in range(0,len(bkgs)+1):                
                    mcErrStrings[bkg] += "\t%.3f"%mcErrs[bkg][i]
            mcErrStrings[bkg]+="\n"
            
    datacard+=binString+processString+processNumberString+rateString+divider
    
    # now nuisances
    datacard+=lumiString
    
    for bkg in bkgs:
        datacard+=mcErrStrings[bkg]

    txtfile = open(txtfileName,"w")
    txtfile.write(datacard)
    txtfile.close()
 
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--lumi-in',dest="lumi_in", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-q','--remove-qcd',dest="removeQCD",default=False,action='store_true',
                  help="remove QCD, while augmenting remaining MC backgrounds")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    lumi = options.lumi
    lumi_in = options.lumi_in
    removeQCD = options.removeQCD
    
    print 'Input files are %s' % ', '.join(args)
    
    #create workspace
    w = rt.RooWorkspace("w"+box)
    variables = initializeWorkspace(w,cfg,box)    
    
    #list of histograms
    ds = []
        
    btagMin =  w.var('nBtag').getMin()
    btagMax =  w.var('nBtag').getMax()

    if removeQCD:
        # first get sum of weights for each background per b-tag bin ( sumW[label] )
        sumW = {}
        sumWQCD = 0.
        for f in args:
            if f.lower().endswith('.root'):
                rootFile = rt.TFile(f)
                tree = rootFile.Get('RazorInclusive')
                if f.lower().find('sms')==-1:
                    
                    label = f.replace('.root','').split('/')[-1]
                    sumW[label] = getSumOfWeights(tree, cfg, box, w, True, f, lumi, lumi_in)
                    if label.find('QCD')!=-1: sumWQCD = sumW[label]
                rootFile.Close()
        # get total sum of weights
        sumWTotal = [sum(allW) for allW in zip( * sumW.values() )]

        # get scale factor to scale other backgrounds by
        k_QCD[box] = [total/(total - qcd) for total, qcd in zip(sumWTotal,sumWQCD)]
         
        print "Sum of Weights Total [ %s ] ="%box, sumWTotal
        print "Sum of Weights QCD   [ %s ] ="%box, sumWQCD
        print "Scale Factor k_QCD   [ %s ] ="%box, k_QCD[box]
    else:        
        z = array('d', cfg.getBinning(box)[2]) # nBtag binning
        k_QCD[box] = [1. for iz in range(1,len(z))]

    #make MC signal and background histograms
    modelString = "" #SMS name
    for i, f in enumerate(args): #loop over input files
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f) #open file
            tree = rootFile.Get('RazorInclusive') #get tree
            if f.lower().find('sms')==-1: #background process
                if removeQCD and f.find('QCD')!=-1:
                    continue # do not add QCD
                else:
                    #set background name according to input file name
                    treeName = ""
                    for name in backgrounds:
                        if f.lower().find(name) != -1:
                            treeName = name
                            break
                    if treeName == "":
                        print("Error: unknown background "+f)
                        sys.exit()
                    #add histogram to output file
                    ds.append(convertTree2TH1(tree, cfg, box, w, True, f, lumi, lumi_in, treeName))
            else: #signal process
                model = f.split('-')[1].split('_')[0]
                massPoint = '_'.join(f.split('_')[3:5])
                modelString = model+'_'+massPoint
                #add histogram to output file
                ds.append(convertTree2TH1(tree, cfg, box, w, True, f ,lumi, lumi_in, modelString))
            rootFile.Close()

    #make data histograms
    #(as a proxy for now, use the sum of the MC)
    data = ds[0].Clone('data_obs')
    for i in range(1, len(ds)): 
        if ds[i].GetName().lower() in backgrounds: 
            data = data + ds[i]
                
    #output file name
    if btagMax>btagMin+1:
        outFileName = 'RazorInclusive_Histograms_lumi-%.1f_%i-%ibtag_%s.root'%(lumi/1000.,btagMin,btagMax-1,box)
    else:
        outFileName = 'RazorInclusive_Histograms_lumi-%.1f_%ibtag_%s.root'%(lumi/1000.,btagMin,box)

    #output histograms
    print "Output File: %s"%(options.outDir+"/"+outFileName)
    outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
    outFile.cd()

    for hist in ds:
        print("Writing histogram: "+hist.GetName())
        hist.Write()
    print("Writing histogram: "+data.GetName())
    data.Write()
   
    outFile.Close()

    #convert dataset list to dict for writing data card
    dsDict = {}
    dsDict["data_obs"] = data
    for d in ds: dsDict[d.GetName()] = d

    #create data card
    writeDataCard_th1(box,modelString,(options.outDir+"/"+outFileName).replace('.root','.txt'),dsDict)
