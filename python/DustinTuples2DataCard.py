from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
import sys
from array import *
from DustinTuple2RooDataSet import initializeWorkspace, getSumOfWeights, boxes, k_T, k_Z, k_W, k_QCD, dPhiCut, MTCut

backgrounds = ['dyjetstoll', 'qcd', 'ttjets', 'zjetstonunu', 'multiboson', 'singletop', 'wjetstolnu', 'ttv']

jet1Cut = 80 #cut on leading jet pt
jet2Cut = 80 #cut on subleading jet pt

def getScaleFactor(tree, treeName, sfs={}, opt=""):
    #get correct value of MR and Rsq
    theMR = tree.MR
    theRsq = tree.Rsq
    if opt == "jesUp": 
        print("Using MR_JESUp")
        theMR = tree.MR_JESUp
        theRsq = tree.Rsq_JESUp
    elif opt == "jesDown": 
        print("Using MR_JESDown")
        theMR = tree.MR_JESDown
        theRsq = tree.Rsq_JESDown
    elif opt == "jerUp": 
        print("Using MR_JERUp")
        theMR = tree.MR_JERUp
        theRsq = tree.Rsq_JERUp
    elif opt == "jerDown": 
        print("Using MR_JERDown")
        theMR = tree.MR_JERDown
        theRsq = tree.Rsq_JERDown

    #get name of scale factor histogram
    centerHistName = ""
    for name in ["ttjets", "wjetstolnu", "dyjetstoll", "zjetstonunu"]:
        if treeName.startswith(name):
            centerHistName = name
            break
    #return 1 if the process has no scale factor histogram
    print("Using centerHistName: "+centerHistName)
    if centerHistName == "": 
        return 1.0

    scaleFactor = sfs[centerHistName].GetBinContent(sfs[centerHistName].FindFixBin(theMR, theRsq))
    print("scaleFactor: "+str(scaleFactor))
    #move up/down by 1 sigma if computing uncertainties
    if "sfstatUp" in opt:
        scaleFactor += sfs[centerHistName].GetBinError(sfs[centerHistName].FindFixBin(theMR, theRsq))
    elif "sfstatDown" in opt:
        scaleFactor -= sfs[centerHistName].GetBinError(sfs[centerHistName].FindFixBin(theMR, theRsq))
    elif "sfsysUp" in opt:
        scaleFactor = sfs[centerHistName+"_SFSysUp"].GetBinContent(sfs[centerHistName+"_SFSysUp"].FindFixBin(theMR, theRsq))
    elif "sfsysDown" in opt:
        scaleFactor = sfs[centerHistName+"_SFSysDown"].GetBinContent(sfs[centerHistName+"_SFSysDown"].FindFixBin(theMR, theRsq))
    print("after possibly shifting: "+str(scaleFactor))

    return scaleFactor
     
def fillRazor3D(tree, hist, weight, btagCutoff, treeName, sfs={}, opt=""):
    """Fill hist for one event, using opt to specify any systematic, etc, that should be applied.
       Returns the weight that was filled."""
    nBTags = min(tree.nBTaggedJets,btagCutoff)

    #multiply weight by appropriate scale factor
    scaleFactor = getScaleFactor(tree, treeName, sfs, opt)
    weight = weight*scaleFactor
                    
    #default
    if opt == "" or opt == "sfstatUp" or opt == "sfstatDown" or opt == "sfsysUp" or opt == "sfsysDown": 
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #muon scale factor up/down
    elif opt == "muoneffUp":
        weight = weight*tree.sf_muonEffUp
        print("weight = "+str(weight))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "muoneffDown":
        weight = weight*tree.sf_muonEffDown
        print("weight = "+str(weight))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #ele scale factor up/down
    elif opt == "eleeffUp":
        weight = weight*tree.sf_eleEffUp
        print("weight = "+str(weight))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "eleeffDown":
        weight = weight*tree.sf_eleEffDown
        print("weight = "+str(weight))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #btag scale factor up/down
    elif opt == "btagUp":
        weight = weight*tree.sf_btagUp
        print("weight = "+str(weight))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "btagDown":
        weight = weight*tree.sf_btagDown
        print("weight = "+str(weight))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #jet energy scale up/down
    elif opt == "jesUp":
        print("jesUp : "+str(tree.MR)+" "+str(tree.MR_JESUp)+" "+str(tree.MR_JESDown))
        hist.Fill(tree.MR_JESUp, tree.Rsq_JESUp, min(tree.nBTaggedJets_JESUp, btagCutoff), weight);
    elif opt == "jesDown":
        print("jesDown : "+str(tree.MR)+" "+str(tree.MR_JESUp)+" "+str(tree.MR_JESDown))
        hist.Fill(tree.MR_JESDown, tree.Rsq_JESDown, min(tree.nBTaggedJets_JESDown, btagCutoff), weight);

    #jet energy resolution up/down
    elif opt == "jerUp":
        print("jerUp : "+str(tree.MR)+" "+str(tree.MR_JERUp)+" "+str(tree.MR_JERDown))
        hist.Fill(tree.MR_JERUp, tree.Rsq_JERUp, min(tree.nBTaggedJets_JERUp, btagCutoff), weight);
    elif opt == "jerDown":
        print("jerDown : "+str(tree.MR)+" "+str(tree.MR_JERUp)+" "+str(tree.MR_JERDown))
        hist.Fill(tree.MR_JERDown, tree.Rsq_JERDown, min(tree.nBTaggedJets_JERDown, btagCutoff), weight);

    else: 
        print("Error in fillRazor3D: option "+opt+" not recognized!")
        sys.exit()

    return weight

def uncorrelate(hists, sysName):
    """Replaces each histogram whose name contains 'sysName' with many copies that represent uncorrelated bin-by-bin systematics"""
    #get all histograms that match the input string
    toUncorrelate = [name for name in hists if sysName in name]
    print("Treating the following distributions as uncorrelated: ")
    for name in toUncorrelate: print name
    
    for name in toUncorrelate:
        print("Uncorrelating "+name)
        #get histogram with central values
        centerName = name.split("_")[:-1]
        centerName = '_'.join(centerName)
        systName = name.split("_")[-1].replace("Up","").replace("Down","")
        print("Central values taken from "+centerName)
        #for each bin create a new histogram in which that bin is up/down and the rest are centered
        for b in range(1,hists[name].GetNbinsX()+1):
            if "Up" in name: 
                newHistName = centerName+"_"+systName+str(b)+"Up"
            elif "Down" in name:
                newHistName = centerName+"_"+systName+str(b)+"Down"
            else: 
                print("Error: shape histogram name "+name+" needs to contain 'Up' or 'Down'")
                return
            hists[newHistName] = hists[centerName].Clone(newHistName)
            hists[newHistName].SetDirectory(0)
            hists[newHistName].SetBinContent(b, hists[name].GetBinContent(b)) #new hist has the unperturbed value in every bin except one
            hists[newHistName].SetBinError(b, hists[name].GetBinError(b))

        #remove the original histogram
        del hists[name]


def convertTree2TH1(tree, cfg, box, workspace, f, lumi, lumi_in, treeName, sfs={}, option=""):
    """Create 1D histogram for direct use with Combine"""
    
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
    myTH3.SetDirectory(0)
    myTH3.Sumw2()

    #temp histogram for btag bins
    htemp = rt.TH1F('htemp_%s'%label,'htemp_%s'%label,len(z)-1,z)

    btagCutoff = 3
    if box in ["MuEle", "MuMu", "EleEle"]:
        btagCutoff = 1
        
    boxCut = boxes[box]
    cuts = 'MR > %f && MR < %f && Rsq > %f && Rsq < %f && min(nBTaggedJets,%i) >= %i && min(nBTaggedJets,%i) < %i && %s && abs(dPhiRazor) < %f && leadingJetPt > %f && subleadingJetPt > %f' % (mRmin,mRmax,rsqMin,rsqMax,btagCutoff,btagMin,btagCutoff,btagMax,boxCut,dPhiCut,jet1Cut,jet2Cut)

    #modify cuts according to box label and/or systematic uncertainty
    if box in ["MuJet", "MuMultiJet", "MuFourJet", "MuSixJet", "EleJet", "EleMultiJet", "EleFourJet", "EleSixJet"]: cuts = cuts+" && mT > "+str(MTCut)
    if box in ["LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet", "LooseLeptonMultiJet"]: cuts = cuts+" && mTLoose > "+str(MTCut)
    if option == "jesUp": 
        cuts = cuts.replace("MR", "MR_JESUp")
        cuts = cuts.replace("Rsq", "Rsq_JESUp")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JESUp")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JESUp")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JESUp")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JESUp")
        cuts = cuts.replace("box", "box_JESUp")
    if option == "jesDown": 
        cuts = cuts.replace("MR", "MR_JESDown")
        cuts = cuts.replace("Rsq", "Rsq_JESDown")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JESDown")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JESDown")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JESDown")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JESDown")
        cuts = cuts.replace("box", "box_JESDown")
    if option == "jerUp": 
        cuts = cuts.replace("MR", "MR_JERUp")
        cuts = cuts.replace("Rsq", "Rsq_JERUp")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JERUp")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JERUp")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JERUp")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JERUp")
        cuts = cuts.replace("box", "box_JERUp")
    if option == "jerDown": 
        cuts = cuts.replace("MR", "MR_JERDown")
        cuts = cuts.replace("Rsq", "Rsq_JERDown")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JERDown")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JERDown")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JERDown")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JERDown")
        cuts = cuts.replace("box", "box_JERDown")

    print("Cuts: "+cuts)

    #get list of entries passing the cuts
    tree.Draw('>>elist', cuts, 'entrylist')
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
        theWeight = tree.weight*lumi*k[btag_bin]/lumi_in
        #if 'SMS' in f: theWeight *= 1
        #else: theWeight *= 100
        filledWeight = fillRazor3D(tree, myTH3, theWeight, btagCutoff, treeName, sfs, option)
        numEntriesByBtag[btag_bin] += 1
        sumEntriesByBtag[btag_bin] += filledWeight

    #unroll into TH1F
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    myTH1 = rt.TH1F(treeName,treeName,nBins,0,nBins)
    myTH1.SetDirectory(0) #prevent it from going out of scope
    myTH1.Sumw2()
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
    lumiErrs = [1.05] #5% lumi systematic on signal
    lumiErrs.extend([1.05 for bkg in bkgs]) #5% lumi systematic on each background
    mcErrs = {} #dictionary of uncorrelated mc bkgd lnN uncertainties

    #get list of shape uncertainties
    shapeNames = []
    for name in hists:
        if "Down" in name:
            shapeName = (name.split("_")[-1]).replace("Down","") #extract the name of the shape histogram
            if shapeName not in shapeNames: shapeNames.append(shapeName)
    #strings listing shape uncertainties for each bkg
    shapeErrs = {name:["1.0"] if model+"_"+name+"Down" in hists else ["-"] for name in shapeNames}
    for name in shapeNames: 
        shapeErrs[name].extend(["1.0" if bkg+"_"+name+"Down" in hists else "-" for bkg in bkgs])

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
    shapeErrStrings = {name:name+"\tshape" for name in shapeNames}
    for name in shapeNames: 
        for i in range(0, len(bkgs)+1):
            shapeErrStrings[name] += "\t"+shapeErrs[name][i]
        shapeErrStrings[name]+="\n"

    datacard+=binString+processString+processNumberString+rateString+divider
    
    # now nuisances
    datacard+=lumiString #lumi uncertainty
    
    for bkg in bkgs:
        datacard+=mcErrStrings[bkg] #MC normalization uncertainties
    for name in shapeNames:
        datacard+=shapeErrStrings[name] #shape uncertainties

    #write card
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
    parser.add_option('--dphi-cut',dest="dPhiCut",default=-1.0,type="float",
                  help="set delta phi cut on the razor hemispheres")
    parser.add_option('--mt-cut',dest="MTCut",default=-1.0,type="float",
                  help="set transverse mass cut")
    parser.add_option('--jet1-cut',dest="jet1Cut",default=-1.0,type="float",
                  help="set leading jet pt cut")
    parser.add_option('--jet2-cut',dest="jet2Cut",default=-1.0,type="float",
                  help="set subleading jet pt cut")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    lumi = options.lumi
    lumi_in = options.lumi_in
    removeQCD = options.removeQCD

    #get appropriate dPhi cut 
    if options.dPhiCut >= 0: dPhiCut = options.dPhiCut
    else:
        if box in ["DiJet", "FourJet", "SixJet", "MultiJet", "MuMu", "MuEle", "EleEle"]:
            dPhiCut = 2.8
        else:
            dPhiCut = 3.2 #no dPhi cut for lepton or loose lepton boxes

    #get appropriate MT cut
    if options.MTCut >= 0: MTCut = options.MTCut
    else:
        if box in ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet", 
                   "EleJet", "EleFourJet", "EleSixJet", "EleMultiJet",
                   "LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet", "LooseLeptonMultiJet"]:
            MTCut = 100 #apply MT > 100 in all lepton and loose lepton boxes
        else:
            MTCut = -1

    if options.jet1Cut >= 0: jet1Cut = options.jet1Cut
    if options.jet2Cut >= 0: jet2Cut = options.jet2Cut

    #get data/MC scale factors from files
    sfFilenames = {
            "ttjets" : "data/ScaleFactors/Placeholders/DummyRun2TTJetsSF.root",
            "wjetstolnu" : "data/ScaleFactors/Placeholders/DummyRun2WJetsSF.root",
            "dyjetstoll" : "data/ScaleFactors/Placeholders/DummyRun2DYJetsSF.root",
            "zjetstonunu" : "data/ScaleFactors/Placeholders/DummyRun2ZNuNuSF.root",
            }
    sfFiles = {name : rt.TFile(sfFilenames[name]) for name in sfFilenames}
    sfHists = {}
    sfHists["ttjets"] = sfFiles["ttjets"].Get("TTJetsSingleLepton")
    sfHists["ttjets_SFSysUp"] = sfFiles["ttjets"].Get("TTJetsSingleLeptonUp")
    sfHists["ttjets_SFSysDown"] = sfFiles["ttjets"].Get("TTJetsSingleLeptonDown")
    sfHists["wjetstolnu"] = sfFiles["wjetstolnu"].Get("WJetsSingleLepton")
    sfHists["wjetstolnu_SFSysUp"] = sfFiles["wjetstolnu"].Get("WJetsSingleLeptonUp")
    sfHists["wjetstolnu_SFSysDown"] = sfFiles["wjetstolnu"].Get("WJetsSingleLeptonDown")
    sfHists["dyjetstoll"] = sfFiles["dyjetstoll"].Get("DYJetsDilepton")
    sfHists["dyjetstoll_SFSysUp"] = sfFiles["dyjetstoll"].Get("DYJetsDileptonUp")
    sfHists["dyjetstoll_SFSysDown"] = sfFiles["dyjetstoll"].Get("DYJetsDileptonDown")
    sfHists["zjetstonunu"] = sfFiles["zjetstonunu"].Get("ZNuNuGJets")
    sfHists["zjetstonunu_SFSysUp"] = sfFiles["zjetstonunu"].Get("ZNuNuGJetsUp")
    sfHists["zjetstonunu_SFSysDown"] = sfFiles["zjetstonunu"].Get("ZNuNuGJetsDown")
    
    #list of shape systematics to apply.
    #if a list of physics processes is given, the uncertainty will be applied to each process in the list, assumed uncorrelated from process to process.
    #if an empty list is given, the uncertainty will be applied (correlated) to all processes, including signal.
    shapes = {
            "muoneff" : [],
            "eleeff" : [],
            "btag" : [],
            "jes" : [],
            "jer" : [],
            "sfstat" : ["ttjets", "wjetstolnu", "dyjetstoll", "zjetstonunu"],
            "sfsys" : ["ttjets", "wjetstolnu", "dyjetstoll", "zjetstonunu"],
            }

    #specify which systematics should be treated as uncorrelated bin-by-bin
    uncorrShapes = [
            #"norm"
            ] 

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
                    print("Building histogram for "+treeName)
                    ds.append(convertTree2TH1(tree, cfg, box, w, f, lumi, lumi_in, treeName, sfs=sfHists))
                    ###get up/down histograms for shape systematics
                    for shape in shapes:
                        for updown in ["Up", "Down"]:
                            if shapes[shape] == []:
                                print("Building histogram for "+treeName+"_"+shape+updown)
                                ds.append(convertTree2TH1(tree, cfg, box, w, f, lumi, lumi_in, treeName+"_"+shape+updown, sfs=sfHists, option=shape+updown))
                            elif treeName.lower() in [s.lower() for s in shapes[shape]]:
                                print("Building histogram for "+treeName+"_"+shape+(treeName.replace('_',''))+updown)
                                ds.append(convertTree2TH1(tree, cfg, box, w, f, lumi, lumi_in, treeName+"_"+shape+(treeName.replace('_',''))+updown, sfs=sfHists, option=shape+updown))
            else: #signal process
                model = f.split('-')[1].split('_')[0]
                massPoint = '_'.join(f.split('_')[3:5])
                modelString = model+'_'+massPoint
                #add histogram to output file
                print("Building histogram for "+modelString)
                ds.append(convertTree2TH1(tree, cfg, box, w, f ,lumi, lumi_in, modelString, sfs=sfHists))
                for shape in shapes:
                    for updown in ["Up", "Down"]:
                        if shapes[shape] == []:
                            print("Building histogram for "+modelString+"_"+shape+updown)
                            ds.append(convertTree2TH1(tree, cfg, box, w, f, lumi, lumi_in, modelString+"_"+shape+updown, sfs=sfHists, option=shape+updown))
                        elif "signal" in [s.lower() for s in shapes[shape]]:
                            print("Building histogram for "+modelString+"_"+shape+(modelString.replace('_',''))+updown)
                            ds.append(convertTree2TH1(tree, cfg, box, w, f, lumi, lumi_in, modelString+"_"+shape+"signal"+updown, sfs=sfHists, option=shape+updown))

            rootFile.Close()

    #convert dataset list to dict
    dsDict = {}
    for d in ds: dsDict[d.GetName()] = d

    #perform uncorrelation procedure 
    for shape in uncorrShapes:
        uncorrelate(dsDict, shape)

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

    #output file
    print "Output File: %s"%(options.outDir+"/"+outFileName)
    outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
    outFile.cd()

    for name in dsDict:
        print("Writing histogram: "+dsDict[name].GetName())
        dsDict[name].Write()
    print("Writing histogram: "+data.GetName())
    data.Write()
   
    outFile.Close()

    #add data for writing card
    dsDict["data_obs"] = data
    #create data card
    writeDataCard_th1(box,modelString,(options.outDir+"/"+outFileName).replace('.root','.txt'),dsDict)
