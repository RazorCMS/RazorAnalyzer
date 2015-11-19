from optparse import OptionParser
import ROOT as rt
import sys
from array import *

#local imports
import rootTools
from framework import Config
from DustinTuple2RooDataSet import initializeWorkspace, getSumOfWeights, boxes, k_T, k_Z, k_W, dPhiCut, MTCut, getCuts
from RunCombine import exec_me

jet1Cut = 80 #cut on leading jet pt
jet2Cut = 80 #cut on subleading jet pt

backgrounds = ['dyjetstoll', 'qcd', 'ttjets', 'zjetstonunu', 'multiboson', 'singletop', 'wjetstolnu', 'ttv']

def getScaleFactor(tree, treeName, sfs={}, opt=""):
    #get correct value of MR and Rsq
    theMR = tree.MR
    theRsq = tree.Rsq
    if opt == "jesUp": 
        theMR = tree.MR_JESUp
        theRsq = tree.Rsq_JESUp
    elif opt == "jesDown": 
        theMR = tree.MR_JESDown
        theRsq = tree.Rsq_JESDown
    elif opt == "jerUp": 
        theMR = tree.MR_JERUp
        theRsq = tree.Rsq_JERUp
    elif opt == "jerDown": 
        theMR = tree.MR_JERDown
        theRsq = tree.Rsq_JERDown
    elif opt == "mesUp":
        theMR = tree.MR_MESUp
        theRsq = tree.Rsq_MESUp
    elif opt == "mesDown":
        theMR = tree.MR_MESDown
        theRsq = tree.Rsq_MESDown
    elif opt == "eesUp":
        theMR = tree.MR_EESUp
        theRsq = tree.Rsq_EESUp
    elif opt == "eesDown":
        theMR = tree.MR_EESDown
        theRsq = tree.Rsq_EESDown

    #get name of scale factor histogram
    centerHistName = ""
    for name in ["ttjets", "wjetstolnu", "dyjetstoll", "zjetstonunu"]:
        if treeName.startswith(name):
            centerHistName = name
            break
    #return 1 if the process has no scale factor histogram
    if centerHistName == "": 
        return 1.0

    scaleFactor = sfs[centerHistName].GetBinContent(sfs[centerHistName].FindFixBin(theMR, theRsq))
    #move up/down by 1 sigma if computing uncertainties
    if "sfstatUp" in opt:
        scaleFactor += sfs[centerHistName].GetBinError(sfs[centerHistName].FindFixBin(theMR, theRsq))
    elif "sfstatDown" in opt:
        scaleFactor -= sfs[centerHistName].GetBinError(sfs[centerHistName].FindFixBin(theMR, theRsq))
    elif "sfsysUp" in opt:
        scaleFactor = sfs[centerHistName+"_sfsysUp"].GetBinContent(sfs[centerHistName+"_sfsysUp"].FindFixBin(theMR, theRsq))
    elif "sfsysDown" in opt:
        scaleFactor = sfs[centerHistName+"_sfsysDown"].GetBinContent(sfs[centerHistName+"_sfsysDown"].FindFixBin(theMR, theRsq))
    elif "sfmethodologyUp" in opt:
        scaleFactor = sfs[centerHistName+"_sfmethodologyUp"].GetBinContent(sfs[centerHistName+"_sfmethodologyUp"].FindFixBin(theMR, theRsq))
    elif "sfmethodologyDown" in opt:
        scaleFactor = sfs[centerHistName+"_sfmethodologyDown"].GetBinContent(sfs[centerHistName+"_sfmethodologyDown"].FindFixBin(theMR, theRsq))

    return scaleFactor
     
def fillRazor3D(tree, hist, weight, btagCutoff, treeName, sfs={}, opt="", sumPdfWeights=None, sumScaleWeights=None, nevents=None):
    """Fill hist for one event, using opt to specify any systematic, etc, that should be applied.
       Returns the weight that was filled."""
    nBTags = min(tree.nBTaggedJets,btagCutoff)

    #check for any required weight histogram
    if 'facscale' in opt or 'renscale' in opt or 'facrenscale' in opt or 'pdf' in opt:
        if nevents is None:
            print "Error in fillRazor3D: no NEvents histogram given!"
            return 
        else:
            integral = nevents.Integral()

    if 'facscale' in opt or 'renscale' in opt or 'facrenscale' in opt:
        if sumScaleWeights is None:
            print "Error in fillRazor3D: no SumScaleWeights histogram given!" 
            return

    if 'pdf' in opt:
        if sumPdfWeights is None:
            print "Error in fillRazor3D: no SumPdfWeights histogram given!"
            return

    #multiply weight by appropriate scale factor
    scaleFactor = getScaleFactor(tree, treeName, sfs, opt)
    weight = weight*scaleFactor
                    
    #default
    if opt == "" or opt == "sfstatUp" or opt == "sfstatDown" or opt == "sfsysUp" or opt == "sfsysDown" or opt == "sfmethodologyUp" or opt == "sfmethodologyDown": 
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #muon scale factor up/down
    elif opt == "muoneffUp":
        weight = weight*tree.sf_muonEffUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "muoneffDown":
        weight = weight*tree.sf_muonEffDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #ele scale factor up/down
    elif opt == "eleeffUp":
        weight = weight*tree.sf_eleEffUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "eleeffDown":
        weight = weight*tree.sf_eleEffDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #muon trig scale factor up/down
    elif opt == "muontrigUp":
        weight = weight*tree.sf_muonTrigUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "muontrigDown":
        weight = weight*tree.sf_muonTrigDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #ele trig scale factor up/down
    elif opt == "eletrigUp":
        weight = weight*tree.sf_eleTrigUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "eletrigDown":
        weight = weight*tree.sf_eleTrigDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #btag scale factor up/down
    elif opt == "btagUp":
        weight = weight*tree.sf_btagUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "btagDown":
        weight = weight*tree.sf_btagDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #muon fastsim scale factor up/down
    elif opt == "muonfastsimUp":
        weight = weight*tree.sf_muonEffFastsimSFUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "muonfastsimDown":
        weight = weight*tree.sf_muonEffFastsimSFDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #ele fastsim scale factor up/down
    elif opt == "elefastsimUp":
        weight = weight*tree.sf_eleEffFastsimSFUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "elefastsimDown":
        weight = weight*tree.sf_eleEffFastsimSFDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #btag fastsim scale factor up/down
    elif opt == "btagfastsimUp":
        weight = weight*tree.sf_btagFastsimSFUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "btagfastsimDown":
        weight = weight*tree.sf_btagFastsimSFDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #ren/fac scale factor up/down
    elif opt == "facscaleUp":
        weight = weight*tree.sf_facScaleUp*integral/sumScaleWeights.GetBinContent(1)
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "facscaleDown":
        weight = weight*tree.sf_facScaleDown*integral/sumScaleWeights.GetBinContent(2)
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "renscaleUp":
        weight = weight*tree.sf_renScaleUp*integral/sumScaleWeights.GetBinContent(3)
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "renscaleDown":
        weight = weight*tree.sf_renScaleDown*integral/sumScaleWeights.GetBinContent(4)
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "facrenscaleUp":
        weight = weight*tree.sf_facRenScaleUp*integral/sumScaleWeights.GetBinContent(5)
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "facrenscaleDown":
        weight = weight*tree.sf_facRenScaleDown*integral/sumScaleWeights.GetBinContent(6)
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #pdf weights
    elif 'pdfUp' in opt:
        pdfNum = int(opt.replace('pdfUp',''))
        weight = weight*(tree.pdfWeights[pdfNum]/tree.genWeight*integral/sumPdfWeights.GetBinContent(pdfNum+1))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif 'pdfDown' in opt:
        pdfNum = int(opt.replace('pdfDown',''))
        weight = weight/(tree.pdfWeights[pdfNum]/tree.genWeight*integral/sumPdfWeights.GetBinContent(pdfNum+1))
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #lumi
    elif 'lumiUp' in opt:
        weight = weight*1.05
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif 'lumiDown' in opt:
        weight = weight/1.05
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #jet energy scale up/down
    elif opt == "jesUp":
        hist.Fill(tree.MR_JESUp, tree.Rsq_JESUp, min(tree.nBTaggedJets_JESUp, btagCutoff), weight)
    elif opt == "jesDown":
        hist.Fill(tree.MR_JESDown, tree.Rsq_JESDown, min(tree.nBTaggedJets_JESDown, btagCutoff), weight)

    #jet energy resolution up/down
    elif opt == "jerUp":
        hist.Fill(tree.MR_JERUp, tree.Rsq_JERUp, min(tree.nBTaggedJets_JERUp, btagCutoff), weight)
    elif opt == "jerDown":
        hist.Fill(tree.MR_JERDown, tree.Rsq_JERDown, min(tree.nBTaggedJets_JERDown, btagCutoff), weight)

    #muon energy scale
    elif opt == 'mesUp':
        hist.Fill(tree.MR_MESUp, tree.Rsq_MESUp, min(tree.nBTaggedJets_MESUp,btagCutoff), weight)
    elif opt == 'mesDown':
        hist.Fill(tree.MR_MESDown, tree.Rsq_MESDown, min(tree.nBTaggedJets_MESDown,btagCutoff), weight)

    #electron energy scale
    elif opt == 'eesUp':
        hist.Fill(tree.MR_EESUp, tree.Rsq_EESUp, min(tree.nBTaggedJets_EESUp,btagCutoff), weight)
    elif opt == 'eesDown':
        hist.Fill(tree.MR_EESDown, tree.Rsq_EESDown, min(tree.nBTaggedJets_EESDown,btagCutoff), weight)


    else: 
        print("Error in fillRazor3D: option "+opt+" not recognized!")
        sys.exit()

    return weight

def uncorrelate(hists, sysName):
    """Replaces each histogram whose name contains 'sysName' with many copies that represent uncorrelated bin-by-bin systematics.
    If referenceHist is given, bins in the shape histograms are assumed to be correlated if and only if they lie in the same bin of referenceHist."""
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

def uncorrelateSFs(hists, sysName, referenceHists, cfg, box):
    """Same as uncorrelate(), but treats bins as correlated if they lie inside the same bin in the reference histogram.
    Needs a config and a box name, to get the correct bin configuration for the razor histogram"""
    #get all histograms that match the input string
    toUncorrelate = [name for name in hists if sysName in name]
    print("Treating the following distributions as uncorrelated: ")
    for name in toUncorrelate: print name

    #make histogram with razor binning
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    myTH3 = rt.TH3D("razor3d","razor3d",len(x)-1,x,len(y)-1,y,len(z)-1,z)
    
    for name in toUncorrelate:
        print("Using reference histogram to determine bin correlations for "+name)
        #get histogram with central values
        centerName = name.split("_")[:-1]
        centerName = '_'.join(centerName)
        systName = name.split("_")[-1].replace("Up","").replace("Down","")
        print("Central values taken from "+centerName)
        #get reference histogram for scale factor binning
        referenceHist = referenceHists[centerName]
        #for each bin create a new histogram in which that bin is up/down and the rest are centered
        for bx in range(1,referenceHist.GetNbinsX()+1):
            for by in range(1,referenceHist.GetNbinsY()+1):
                b = referenceHist.GetBin(bx,by)
                if "Up" in name: 
                    newHistName = centerName+"_"+systName+str(b)+"Up"
                elif "Down" in name:
                    newHistName = centerName+"_"+systName+str(b)+"Down"
                else: 
                    print("Error: shape histogram name "+name+" needs to contain 'Up' or 'Down'")
                    return
                hists[newHistName] = hists[centerName].Clone(newHistName)
                hists[newHistName].SetDirectory(0)
                #find bins in hists[name] that lie inside bin b of referenceHist
                i = 0
                for ix in range(1,len(x)):
                    for iy in range(1,len(y)):
                        for iz in range(1,len(z)):
                            #i = 1D histogram bin index
                            i+= 1
                            #get MR and Rsq at center of bin in 3d histogram
                            mrCenter = myTH3.GetXaxis().GetBinCenter(ix)
                            rsqCenter = myTH3.GetYaxis().GetBinCenter(iy)
                            #correct MR or Rsq if they lie outside the reference histogram
                            if mrCenter > referenceHist.GetXaxis().GetXmax(): 
                                mrCenter = referenceHist.GetXaxis().GetXmax() - 1
                            if rsqCenter > referenceHist.GetYaxis().GetXmax():
                                rsqCenter = referenceHist.GetYaxis().GetXmax() - 0.01
                            #if the bin matches the current reference histogram bin, update the contents
                            if referenceHist.FindFixBin(mrCenter, rsqCenter) == b: #bin matches
                                hists[newHistName].SetBinContent(i, hists[name].GetBinContent(i)) #new hist has the unperturbed value in every bin except one
                                hists[newHistName].SetBinError(i, hists[name].GetBinError(i))

        #remove the original histogram
        del hists[name]

def convertTree2TH1(tree, cfg, box, workspace, f, globalScaleFactor, treeName, sfs={}, sysErrOpt="", pileupWeightHist=None, hadronicTriggerWeight=None, sumPdfWeights=None, sumScaleWeights=None, nevents=None):
    """Create 1D histogram for direct use with Combine"""
    
    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    #get k factor 
    if 'TTJets' in f:
        k = k_T
    elif 'DYJets' in f or 'ZJets' in f:
        k = k_Z
    elif 'WJets' in f:
        k = k_W
    else:
        k = 1.
    btagCutoff = 3

    #make histogram with razor binning
    label = f.replace('.root','').split('/')[-1]
    myTH3 = rt.TH3D(treeName+"3d",treeName+"3d",len(x)-1,x,len(y)-1,y,len(z)-1,z)
    myTH3.SetDirectory(0)
    myTH3.Sumw2()

    #temp histogram for btag bins
    htemp = rt.TH1F('htemp_%s'%label,'htemp_%s'%label,len(z)-1,z)

    cuts = getCuts(workspace, box)

    #modify cuts based on histogram option
    if sysErrOpt == "jesUp": 
        cuts = cuts.replace("MR", "MR_JESUp")
        cuts = cuts.replace("Rsq", "Rsq_JESUp")
        cuts = cuts.replace("mT ", "mT_JESUp ")
        cuts = cuts.replace("mTLoose", "mTLoose_JESUp")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JESUp")
        cuts = cuts.replace("nJets80", "nJets80_JESUp")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JESUp")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JESUp")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JESUp")
        cuts = cuts.replace("box", "box_JESUp")
    if sysErrOpt == "jesDown": 
        cuts = cuts.replace("MR", "MR_JESDown")
        cuts = cuts.replace("Rsq", "Rsq_JESDown")
        cuts = cuts.replace("mT ", "mT_JESDown ")
        cuts = cuts.replace("mTLoose", "mTLoose_JESDown")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JESDown")
        cuts = cuts.replace("nJets80", "nJets80_JESDown")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JESDown")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JESDown")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JESDown")
        cuts = cuts.replace("box", "box_JESDown")
    if sysErrOpt == "jerUp": 
        cuts = cuts.replace("MR", "MR_JERUp")
        cuts = cuts.replace("Rsq", "Rsq_JERUp")
        cuts = cuts.replace("mT ", "mT_JERUp ")
        cuts = cuts.replace("mTLoose", "mTLoose_JERUp")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JERUp")
        cuts = cuts.replace("nJets80", "nJets80_JERUp")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JERUp")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JERUp")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JERUp")
        cuts = cuts.replace("box", "box_JERUp")
    if sysErrOpt == "jerDown": 
        cuts = cuts.replace("MR", "MR_JERDown")
        cuts = cuts.replace("Rsq", "Rsq_JERDown")
        cuts = cuts.replace("mT ", "mT_JERUp ")
        cuts = cuts.replace("mTLoose", "mTLoose_JERUp")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_JERDown")
        cuts = cuts.replace("nJets80", "nJets80_JERDown")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_JERDown")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_JERDown")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_JERDown")
        cuts = cuts.replace("box", "box_JERDown")
    if sysErrOpt == "mesUp": 
        cuts = cuts.replace("MR", "MR_MESUp")
        cuts = cuts.replace("Rsq", "Rsq_MESUp")
        cuts = cuts.replace("mT ", "mT_MESUp ")
        cuts = cuts.replace("mTLoose", "mTLoose_MESUp")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_MESUp")
        cuts = cuts.replace("nJets80", "nJets80_MESUp")
        cuts = cuts.replace('leadingTightMuPt', 'leadingTightMuPt_MESUp')
        cuts = cuts.replace("nSelectedJets", "nSelectedJets_MESUp")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_MESUp")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_MESUp")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_MESUp")
        cuts = cuts.replace("box", "box_MESUp")
    if sysErrOpt == "mesDown": 
        cuts = cuts.replace("MR", "MR_MESDown")
        cuts = cuts.replace("Rsq", "Rsq_MESDown")
        cuts = cuts.replace("mT ", "mT_MESDown ")
        cuts = cuts.replace("mTLoose", "mTLoose_MESDown")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_MESDown")
        cuts = cuts.replace("nJets80", "nJets80_MESDown")
        cuts = cuts.replace('leadingTightMuPt', 'leadingTightMuPt_MESDown')
        cuts = cuts.replace("nSelectedJets", "nSelectedJets_MESDown")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_MESDown")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_MESDown")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_MESDown")
        cuts = cuts.replace("box", "box_MESDown")
    if sysErrOpt == "eesUp": 
        cuts = cuts.replace("MR", "MR_EESUp")
        cuts = cuts.replace("Rsq", "Rsq_EESUp")
        cuts = cuts.replace("mT ", "mT_EESUp ")
        cuts = cuts.replace("mTLoose", "mTLoose_EESUp")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_EESUp")
        cuts = cuts.replace("nJets80", "nJets80_EESUp")
        cuts = cuts.replace('leadingTightMuPt', 'leadingTightMuPt_EESUp')
        cuts = cuts.replace("nSelectedJets", "nSelectedJets_EESUp")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_EESUp")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_EESUp")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_EESUp")
        cuts = cuts.replace("box", "box_EESUp")
    if sysErrOpt == "eesDown": 
        cuts = cuts.replace("MR", "MR_EESDown")
        cuts = cuts.replace("Rsq", "Rsq_EESDown")
        cuts = cuts.replace("mT ", "mT_EESDown ")
        cuts = cuts.replace("mTLoose", "mTLoose_EESDown")
        cuts = cuts.replace("dPhiRazor", "dPhiRazor_EESDown")
        cuts = cuts.replace("nJets80", "nJets80_EESDown")
        cuts = cuts.replace('leadingTightMuPt', 'leadingTightMuPt_EESDown')
        cuts = cuts.replace("nSelectedJets", "nSelectedJets_EESDown")
        cuts = cuts.replace("nBTaggedJets", "nBTaggedJets_EESDown")
        cuts = cuts.replace(" leadingJetPt", " leadingJetPt_EESDown")
        cuts = cuts.replace(" subleadingJetPt", " subleadingJetPt_EESDown")
        cuts = cuts.replace("box", "box_EESDown")

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
        theWeight = tree.weight*k*globalScaleFactor
        #########################
        #temporary reweighting for pileup and hadronic trigger
        if pileupWeightHist is not None:
            pileupWeight = pileupWeightHist.GetBinContent(pileupWeightHist.GetXaxis().FindFixBin(tree.nVtx))
            theWeight *= pileupWeight
        if hadronicTriggerWeight is not None:
            if box in ["MultiJet", "DiJet", "FourJet", "SixJet", "LooseLeptonMultiJet", "LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet"]:
                theWeight *= hadronicTriggerWeight
        #########################
        filledWeight = fillRazor3D(tree, myTH3, theWeight, btagCutoff, treeName, sfs, sysErrOpt, sumPdfWeights=sumPdfWeights, sumScaleWeights=sumScaleWeights, nevents=nevents)
        numEntriesByBtag[btag_bin] += 1
        sumEntriesByBtag[btag_bin] += filledWeight

    #unroll into TH1F
    nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
    maxBins = 244
    myTH1 = rt.TH1F(treeName,treeName,maxBins,0,maxBins)
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

    #10% normalization uncertainty on signal (proxy for cross section error)
    mcErrs[model] = [1.10]
    for bkg in bkgs:
        mcErrs[model].extend([1.00 for bkg1 in bkgs])
    #10% normalization uncertainty on rare backgrounds
    for bkg in bkgs:
        mcErrs[bkg] = [1.00]
        mcErrs[bkg].extend([1.00 + 0.10*(bkg==bkg1 and bkg1 not in 
                ['ttjets','wjetstolnu','dyjetstoll','zjetstonunu']) for bkg1 in bkgs]) 
            
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
    mcErrStrings[model] = "%s_norm\tlnN"%(model)
    for i in range(0,len(bkgs)+1):                
        mcErrStrings[model] += "\t%.3f"%mcErrs[model][i]
    mcErrStrings[model]+="\n"
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
    
    datacard+=mcErrStrings[model]
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
    boxList = box.split('_')
    lumi = options.lumi
    lumi_in = options.lumi_in

    for curBox in boxList:
        #get appropriate dPhi cut 
        if options.dPhiCut >= 0: dPhiCut = options.dPhiCut
        else:
            if curBox in ["DiJet", "FourJet", "SixJet", "MultiJet", "MuMu", "MuEle", "EleEle"]:
                dPhiCut = 2.8
            else:
                dPhiCut = 3.2 #no dPhi cut for lepton or loose lepton boxes

        #get appropriate MT cut
        if options.MTCut >= 0: MTCut = options.MTCut
        else:
            if curBox in ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet", 
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
        sfHists["ttjets_sfsysUp"] = sfFiles["ttjets"].Get("TTJetsSingleLeptonUp")
        sfHists["ttjets_sfsysDown"] = sfFiles["ttjets"].Get("TTJetsSingleLeptonDown")
        sfHists["wjetstolnu"] = sfFiles["wjetstolnu"].Get("WJetsSingleLepton")
        sfHists["wjetstolnu_sfsysUp"] = sfFiles["wjetstolnu"].Get("WJetsSingleLeptonUp")
        sfHists["wjetstolnu_sfsysDown"] = sfFiles["wjetstolnu"].Get("WJetsSingleLeptonDown")
        sfHists["dyjetstoll"] = sfFiles["dyjetstoll"].Get("DYJetsDilepton")
        sfHists["dyjetstoll_sfsysUp"] = sfFiles["dyjetstoll"].Get("DYJetsDileptonUp")
        sfHists["dyjetstoll_sfsysDown"] = sfFiles["dyjetstoll"].Get("DYJetsDileptonDown")
        sfHists["zjetstonunu"] = sfFiles["zjetstonunu"].Get("ZNuNuGJets")
        sfHists["zjetstonunu_sfsysUp"] = sfFiles["zjetstonunu"].Get("ZNuNuGJetsUp")
        sfHists["zjetstonunu_sfsysDown"] = sfFiles["zjetstonunu"].Get("ZNuNuGJetsDown")
        #up histogram for the methodology check is ZNuNuDilepton
        sfHists["zjetstonunu_sfmethodologyUp"] = sfFiles["zjetstonunu"].Get("ZNuNuDilepton")
        #down histogram for the methodology check is ZNuNuGJets - (ZNuNuDilepton - ZNuNuGJets)
        sfHists["zjetstonunu_sfmethodologyDown"] = sfHists["zjetstonunu"].Clone("ZNuNuSFMethodDown")
        sfHists["zjetstonunu_sfmethodologyDown"].Add(sfHists["zjetstonunu"])
        sfHists["zjetstonunu_sfmethodologyDown"].Add(sfHists["zjetstonunu_sfmethodologyUp"], -1)

        #list of shape systematics to apply.
        #if a list of physics processes is given, the uncertainty will be applied to each process in the list, assumed uncorrelated from process to process.
        #if an empty list is given, the uncertainty will be applied (correlated) to all processes, including signal.
        shapes = {
                #"muoneff" : [],
                #"eleeff" : [],
                #"btag" : [],
                #"jes" : [],
                #"jer" : [],
                #"sfstat" : ["ttjets", "wjetstolnu", "dyjetstoll", "zjetstonunu"],
                #"sfsys" : ["ttjets", "wjetstolnu", "dyjetstoll", "zjetstonunu"],
                #"sfmethodology" : ["zjetstonunu"]
                }

        #specify which systematics should be treated as uncorrelated bin-by-bin
        uncorrSFShapes = [
                #"sfstat",
                ] 

        print 'Input files are %s' % ', '.join(args)
        
        #create workspace
        w = rt.RooWorkspace("w"+curBox)
        variables = initializeWorkspace(w,cfg,curBox)    
        
        #list of histograms
        ds = []
            
        btagMin =  w.var('nBtag').getMin()
        btagMax =  w.var('nBtag').getMax()
        z = array('d', cfg.getBinning(curBox)[2]) # nBtag binning

        #make MC signal and background histograms
        modelString = "" #SMS name
        for i, f in enumerate(args): #loop over input files
            if f.lower().endswith('.root'):
                rootFile = rt.TFile(f) #open file
                tree = rootFile.Get('RazorInclusive') #get tree
                if f.lower().find('sms')==-1: #background process
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
                    ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=treeName, sfs=sfHists))
                    ###get up/down histograms for shape systematics
                    for shape in shapes:
                        for updown in ["Up", "Down"]:
                            if shapes[shape] == []:
                                print("Building histogram for "+treeName+"_"+shape+updown)
                                ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=treeName+"_"+shape+updown, sfs=sfHists, sysErrOpt=shape+updown))
                            elif treeName.lower() in [s.lower() for s in shapes[shape]]:
                                print("Building histogram for "+treeName+"_"+shape+(treeName.replace('_',''))+updown)
                                ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=treeName+"_"+shape+(treeName.replace('_',''))+updown, sfs=sfHists, sysErrOpt=shape+updown))
                else: #signal process
                    model = f.split('-')[1].split('_')[0]
                    massPoint = '_'.join(f.split('_')[3:5])
                    modelString = model+'_'+massPoint
                    #add histogram to output file
                    print("Building histogram for "+modelString)
                    ds.append(convertTree2TH1(tree, cfg, curBox, w, f , globalScaleFactor=lumi/lumi_in, treeName=modelString, sfs=sfHists))
                    for shape in shapes:
                        for updown in ["Up", "Down"]:
                            if shapes[shape] == []:
                                print("Building histogram for "+modelString+"_"+shape+updown)
                                ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=modelString+"_"+shape+updown, sfs=sfHists, sysErrOpt=shape+updown))
                            elif "signal" in [s.lower() for s in shapes[shape]]:
                                print("Building histogram for "+modelString+"_"+shape+(modelString.replace('_',''))+updown)
                                ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=modelString+"_"+shape+"signal"+updown, sfs=sfHists, sysErrOpt=shape+updown))

                rootFile.Close()

        #convert dataset list to dict
        dsDict = {}
        for d in ds: dsDict[d.GetName()] = d

        #perform uncorrelation procedure
        for shape in uncorrSFShapes:
            uncorrelateSFs(dsDict, shape, sfHists, cfg, curBox)

        #make data histograms
        #(as a proxy for now, use the sum of the MC)
        data = ds[0].Clone('data_obs')
        for i in range(1, len(ds)): 
            if ds[i].GetName().lower() in backgrounds: 
                data = data + ds[i]

        #output file name
        if btagMax>btagMin+1:
            outFileName = 'RazorInclusive_Histograms_lumi-%.1f_%i-%ibtag_%s.root'%(lumi/1000.,btagMin,btagMax-1,curBox)
        else:
            outFileName = 'RazorInclusive_Histograms_lumi-%.1f_%ibtag_%s.root'%(lumi/1000.,btagMin,curBox)

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
        writeDataCard_th1(curBox,modelString,(options.outDir+"/"+outFileName).replace('.root','.txt'),dsDict)

    #run combine
    if len(boxList) == 1:
        #get card name
        if btagMax>btagMin+1:
            cardName = '%s/RazorInclusive_Histograms_lumi-%.1f_%i-%ibtag_%s.txt'%(options.outDir,lumi/1000.,btagMin,btagMax-1,boxList[0])
        else:
            cardName = '%s/RazorInclusive_Histograms_lumi-%.1f_%ibtag_%s.txt'%(options.outDir,lumi/1000.,btagMin,boxList[0])

        exec_me('combine -M Asymptotic '+cardName, False)

    elif len(boxList) > 1:
        #get card names
        cardNames = []
        for curBox in boxList:
            if btagMax>btagMin+1:
                cardName = 'RazorInclusive_Histograms_lumi-%.1f_%i-%ibtag_%s.txt'%(lumi/1000.,btagMin,btagMax-1,curBox)
            else:
                cardName = 'RazorInclusive_Histograms_lumi-%.1f_%ibtag_%s.txt'%(lumi/1000.,btagMin,curBox)
            cardNames.append(cardName)
        #combine cards
        exec_me('cd '+options.outDir+'; combineCards.py '+(' '.join(cardNames))+' > RazorInclusive_Histograms_'+('_'.join(boxList))+'.txt; cd ..', False)
        #call combine
        exec_me('combine -M Asymptotic '+options.outDir+'/RazorInclusive_Histograms_'+('_'.join(boxList))+'.txt', False)


