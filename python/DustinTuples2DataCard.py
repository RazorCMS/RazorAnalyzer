import os
import argparse
import ROOT as rt
import sys
from array import *

#local imports
import rootTools
from framework import Config
from DustinTuple2RooDataSet import initializeWorkspace, boxes, k_T, k_Z, k_W, getCuts
from RunCombine import exec_me
from macro.razorWeights import loadScaleFactorHists
import macro.macro as macro
from MixedBranchingRatioWeight import getBRWeight

DIR_MC = "SimpleBackgrounds"
#DIR_MC = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106/forfit"
DIR_DATA = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p23_ForPreappFreezing20151106"
DIR_SIGNAL = "root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p24_ForApproval20151208/jobs/combined/"
DATA_NAMES={
    'MultiJet':DIR_DATA+'/RazorInclusive_HTMHT_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered.root',
    'EleMultiJet':DIR_DATA+'/RazorInclusive_SingleElectron_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered.root',
    'MuMultiJet':DIR_DATA+'/RazorInclusive_SingleMuon_Run2015D_2093pb_GoodLumiGolden_RazorSkim_Filtered.root',
    }
FILENAMES_MC = {
        "ttjets"        : DIR_MC+"/"+"RazorInclusive_TTJets_Madgraph_1pb_weighted.root",
        "wjetstolnu"    : DIR_MC+"/"+"RazorInclusive_WJetsToLNu_1pb_weighted.root",
        "singletop"     : DIR_MC+"/"+"RazorInclusive_ST_1pb_weighted.root",
        "other"         : DIR_MC+"/"+"RazorInclusive_Other_1pb_weighted.root",
        "dyjetstoll"    : DIR_MC+"/"+"RazorInclusive_DYJetsToLL_1pb_weighted.root",
        "zjetstonunu"   : DIR_MC+"/"+"RazorInclusive_ZJetsToNuNu_1pb_weighted.root",
        "qcd"           : DATA_NAMES["MultiJet"]
        }
#FILENAMES_MC = {
#        "ttjets"       : DIR_MC+"/"+"RazorInclusive_TTJets_Madgraph_Leptonic_1pb_weighted_RazorSkim.root",
#        "wjetstolnu"   : DIR_MC+"/"+"RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted_RazorSkim.root",
#        "singletop"    : DIR_MC+"/"+"RazorInclusive_ST_1pb_weighted_RazorSkim.root",
#        "other"        : DIR_MC+"/"+"RazorInclusive_Other_1pb_weighted_RazorSkim.root",
#        "dyjetstoll"   : DIR_MC+"/"+"RazorInclusive_DYJetsToLL_M-5toInf_HTBinned_1pb_weighted_RazorSkim.root",
#        "zjetstonunu"  : DIR_MC+"/"+"RazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted_RazorSkim.root",
#         "qcd"         : DATA_NAMES["MultiJet"]
#        }

lumiUncertainty = 0.027

def getTheoryCrossSectionAndError(mGluino=-1, mStop=-1):
    thyXsec = -1
    thyXsecErr = -1

    if mGluino!=-1:
        for line in open('data/gluino13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mGluino))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1]) #pb
                thyXsecErr = 0.01*float(line.split(',')[2])
    if mStop!=-1:
        for line in open('data/stop13TeV.txt','r'):
            line = line.replace('\n','')
            if str(int(mStop))==line.split(',')[0]:
                thyXsec = float(line.split(',')[1]) #pb
                thyXsecErr = 0.01*float(line.split(',')[2]) 

    return thyXsec,thyXsecErr

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
    for name in ["ttjets", "wjetstolnu", "zjetstonunu"]:
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
        scaleFactor = sfs[centerHistName+"Up"].GetBinContent(sfs[centerHistName+"Up"].FindFixBin(theMR, theRsq))
    elif "sfsysDown" in opt:
        scaleFactor = sfs[centerHistName+"Down"].GetBinContent(sfs[centerHistName+"Down"].FindFixBin(theMR, theRsq))
    #ttjets cross check: apply extra TTJets scale factor
    elif "ttcrosscheckUp" in opt:
        scaleFactor *= sfs["ttjetsdilepton"].GetBinContent(sfs["ttjetsdilepton"].FindFixBin(theMR, theRsq))
    elif "ttcrosscheckDown" in opt:
        scaleFactor /= sfs["ttjetsdilepton"].GetBinContent(sfs["ttjetsdilepton"].FindFixBin(theMR, theRsq))
    #veto lepton cross check
    elif opt == "vetolepcrosscheckUp":
        if tree.leadingGenLeptonPt > 0 and abs(tree.leadingGenLeptonType) != 15:
            scaleFactor *= sfs["vetolepton"].GetBinContent(sfs["vetolepton"].FindFixBin(tree.leadingGenLeptonPt))
    elif opt == "vetolepcrosscheckDown":
        if tree.leadingGenLeptonPt > 0 and abs(tree.leadingGenLeptonType) != 15:
            scaleFactor /= sfs["vetolepton"].GetBinContent(sfs["vetolepton"].FindFixBin(tree.leadingGenLeptonPt))
    #veto tau cross check
    elif opt == "vetotaucrosscheckUp":
        if tree.leadingGenLeptonPt > 0 and abs(tree.leadingGenLeptonType) == 15:
            scaleFactor *= sfs["vetotau"].GetBinContent(sfs["vetotau"].FindFixBin(tree.leadingGenLeptonPt))
    elif opt == "vetotaucrosscheckDown":
        if tree.leadingGenLeptonPt > 0 and abs(tree.leadingGenLeptonType) == 15:
            scaleFactor /= sfs["vetotau"].GetBinContent(sfs["vetotau"].FindFixBin(tree.leadingGenLeptonPt))

    return scaleFactor
     
def fillRazor3D(tree, hist, weight, btagCutoff, treeName, sfs={}, opt="", sumPdfWeights=None, sumScaleWeights=None, nevents=None):
    """Fill hist for one event, using opt to specify any systematic, etc, that should be applied.
       Returns the weight that was filled."""
    nBTags = min(tree.nBTaggedJets,btagCutoff)

    #check for any required weight histogram
    if 'facscale' in opt or 'renscale' in opt or 'facrenscale' in opt:
    #if 'facscale' in opt or 'renscale' in opt or 'facrenscale' in opt or 'pdf' in opt:
        if nevents is None:
            print "Error in fillRazor3D: no NEvents histogram given!"
            return 
        else:
            integral = nevents.Integral()

    if 'facscale' in opt or 'renscale' in opt or 'facrenscale' in opt:
        if sumScaleWeights is None:
            print "Error in fillRazor3D: no SumScaleWeights histogram given!" 
            return

    #if 'pdf' in opt:
    #    if sumPdfWeights is None:
    #        print "Error in fillRazor3D: no SumPdfWeights histogram given!"
    #        return

    #multiply weight by appropriate scale factor
    scaleFactor = getScaleFactor(tree, treeName, sfs, opt)
    weight = weight*scaleFactor
                    
    #default
    if opt == "" or opt == "sfstatUp" or opt == "sfstatDown" or opt == "sfsysUp" or opt == "sfsysDown" or "mcstat" in opt or "pdf" in opt or 'ttcrosscheck' in opt or "vetolepcrosscheck" in opt or "vetotaucrosscheck" in opt: 
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #muon scale factor up/down
    elif opt == "tightmuoneffUp":
        weight = weight*tree.sf_muonEffUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "tightmuoneffDown":
        weight = weight*tree.sf_muonEffDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "vetomuoneffUp":
        weight = weight*tree.sf_vetoMuonEffUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "vetomuoneffDown":
        weight = weight*tree.sf_vetoMuonEffDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #ele scale factor up/down
    elif opt == "tighteleeffUp":
        weight = weight*tree.sf_eleEffUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "tighteleeffDown":
        weight = weight*tree.sf_eleEffDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "vetoeleeffUp":
        weight = weight*tree.sf_vetoEleEffUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "vetoeleeffDown":
        weight = weight*tree.sf_vetoEleEffDown
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
    elif opt == "mistagUp":
        weight = weight*tree.sf_mistagUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "mistagDown":
        weight = weight*tree.sf_mistagDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #tight muon fastsim scale factor up/down
    elif opt == "tightmuonfastsimUp":
        weight = weight*tree.sf_muonEffFastsimSFUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "tightmuonfastsimDown":
        weight = weight*tree.sf_muonEffFastsimSFDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #tight ele fastsim scale factor up/down
    elif opt == "tightelefastsimUp":
        weight = weight*tree.sf_eleEffFastsimSFUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "tightelefastsimDown":
        weight = weight*tree.sf_eleEffFastsimSFDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #veto muon fastsim scale factor up/down
    elif opt == "vetomuonfastsimUp":
        weight = weight*tree.sf_vetoMuonEffFastsimSFUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == "vetomuonfastsimDown":
        weight = weight*tree.sf_vetoMuonEffFastsimSFDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #veto ele fastsim scale factor up/down
    elif opt == "vetoelefastsimUp":
        weight = weight*tree.sf_vetoEleEffFastsimSFUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight) 
    elif opt == "vetoelefastsimDown":
        weight = weight*tree.sf_vetoEleEffFastsimSFDown
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
    #elif 'pdfUp' in opt:
    #    pdfNum = int(opt.replace('pdfUp','').replace('n',''))
    #    weight = weight*(tree.pdfWeights[pdfNum]/tree.genWeight*integral/sumPdfWeights.GetBinContent(pdfNum+1))
    #    hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    #elif 'pdfDown' in opt:
    #    pdfNum = int(opt.replace('pdfDown','').replace('n',''))
    #    weight = weight/(tree.pdfWeights[pdfNum]/tree.genWeight*integral/sumPdfWeights.GetBinContent(pdfNum+1))
    #    hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #lumi
    elif 'lumiUp' in opt:
        weight = weight*(1+lumiUncertainty)
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif 'lumiDown' in opt:
        weight = weight/(1+lumiUncertainty)
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

    #pileup 
    elif opt == 'pileupUp':
        weight = weight*tree.pileupWeightUp 
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == 'pileupDown':
        weight = weight*tree.pileupWeightDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    #ISR 
    elif opt == 'isrUp':
        weight = weight*tree.ISRSystWeightUp
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)
    elif opt == 'isrDown':
        weight = weight*tree.ISRSystWeightDown
        hist.Fill(tree.MR, tree.Rsq, nBTags, weight)

    else: 
        print("Error in fillRazor3D: option "+opt+" not recognized!")
        sys.exit()

    return weight

def uncorrelate(hists, sysName, suppressLevel=None):
    """Replaces each histogram whose name contains 'sysName' with many copies that represent uncorrelated bin-by-bin systematics.
    suppressLevel: if provided, new histograms will only be created for bins that differ from nominal by a fractional amount greater than suppressLevel."""
    #get all histograms that match the input string
    toUncorrelate = [name for name in hists if sysName in name]
    print "Treating the following distributions as uncorrelated for",sysName,": "
    for name in toUncorrelate: print name
    
    #get names of individual systematics
    systNames = []
    for name in toUncorrelate:
        systName = name.replace("Up","").replace("Down","")
        if systName not in systNames:
            systNames.append(systName)

    for name in systNames:
        print("Uncorrelating "+name)
        #get histogram with central values
        centerName = name.split("_")[:-1]
        centerName = '_'.join(centerName)
        #get up and down variants
        upName = name+'Up'
        downName = name+'Down'
        uncName = name.split("_")[-1]
        print("Central values taken from "+centerName)
        #for each bin create a new histogram in which that bin is up/down and the rest are centered
        for b in range(1,hists[centerName].GetNbinsX()+1):
            newHistUpName = centerName+"_"+uncName+str(b)+"Up"
            newHistDownName = centerName+"_"+uncName+str(b)+"Down"

            #check level of agreement with the nominal
            if suppressLevel is not None:
                #get percent difference from nominal
                if hists[centerName].GetBinContent(b) > 0:
                    percDifferenceUp = abs(hists[upName].GetBinContent(b)-hists[centerName].GetBinContent(b))/hists[centerName].GetBinContent(b)
                    percDifferenceDown = abs(hists[downName].GetBinContent(b)-hists[centerName].GetBinContent(b))/hists[centerName].GetBinContent(b)
                    percDifference = max(percDifferenceUp, percDifferenceDown)
                    if percDifference <= suppressLevel: 
                        #print "Suppressing nuisance in bin",b,"(agrees at",percDifference,"level)"
                        continue
                elif hists[upName].GetBinContent(b) == hists[centerName].GetBinContent(b) and hists[downName].GetBinContent(b) == hists[centerName].GetBinContent(b): 
                        #print "Suppressing nuisance in bin",b,"because there is no change from the nominal"
                        continue

            #new up histogram
            hists[newHistUpName] = hists[centerName].Clone(newHistUpName)
            hists[newHistUpName].SetDirectory(0)
            hists[newHistUpName].SetBinContent(b, hists[upName].GetBinContent(b)) #new hist has the unperturbed value in every bin except one
            hists[newHistUpName].SetBinError(b, hists[upName].GetBinError(b))

            #new down histogram
            hists[newHistDownName] = hists[centerName].Clone(newHistDownName)
            hists[newHistDownName].SetDirectory(0)
            hists[newHistDownName].SetBinContent(b, hists[downName].GetBinContent(b)) #new hist has the unperturbed value in every bin except one
            hists[newHistDownName].SetBinError(b, hists[downName].GetBinError(b))

        #remove the original histogram
        del hists[upName]
        del hists[downName]

def makeNewHistogramForUncorrelateSFs(name, centerName, systName, number, hists):
    """Set up a new histogram, for use with the uncorrelateSFs function"""
    if "Up" in name: 
        newHistName = centerName+"_"+systName+str(number)+"Up"
    elif "Down" in name:
        newHistName = centerName+"_"+systName+str(number)+"Down"
    else: 
        print("Error: shape histogram name "+name+" needs to contain 'Up' or 'Down'")
        return
    hists[newHistName] = hists[centerName].Clone(newHistName)
    hists[newHistName].SetDirectory(0)
    return newHistName

def setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN, sigBN, sysHist, newHist, referenceHist):
    """If the signal bin is inside the reference bin, perturb the signal bin"""
    #correct MR or Rsq if they lie outside the reference histogram
    if mrCenter > referenceHist.GetXaxis().GetXmax(): 
        mrCenter = referenceHist.GetXaxis().GetXmax() - 1
    if rsqCenter > referenceHist.GetYaxis().GetXmax():
        rsqCenter = referenceHist.GetYaxis().GetXmax() - 0.01
    #if the bin matches the current reference histogram bin, update the contents
    if referenceHist.FindBin(mrCenter, rsqCenter) == refBN: #bin matches
        newHist.SetBinContent(sigBN, sysHist.GetBinContent(sigBN)) #new hist has the unperturbed value in every bin except one
        newHist.SetBinError(sigBN, sysHist.GetBinError(sigBN))
        return True
    else:
        return False

def uncorrelateSFs(hists, sysName, referenceHists, cfg, box, unrollBins=None):
    """Same as uncorrelate(), but treats bins as correlated if they lie inside the same bin in the reference histogram.
    Needs a config and a box name, to get the correct bin configuration for the razor histogram"""
    #get all histograms that match the input string
    toUncorrelate = [name for name in hists if sysName in name]
    print "Uncorrelate SFs:",sysName
    print("Treating the following distributions as uncorrelated: ")
    for name in toUncorrelate: print name

    x = array('d', cfg.getBinning(box)[0]) # MR binning
    y = array('d', cfg.getBinning(box)[1]) # Rsq binning
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    #make histogram with razor binning
    myTH3 = rt.TH3D("razor3d"+name,"razor3d",len(x)-1,x,len(y)-1,y,len(z)-1,z)

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
        if referenceHist.InheritsFrom("TH2Poly"):
            for bn in range(1,referenceHist.GetNumberOfBins()+1):
                #print "In bin",bn,"of scale factor histogram"
                matchedAtLeastOneBin = False
                newHistName = makeNewHistogramForUncorrelateSFs(name, centerName, systName, bn, hists)
                #find all bins of signal histogram that are within this bin
                if unrollBins is None:
                    i = 0
                    for ix in range(1,len(x)):
                        for iy in range(1,len(y)):
                            for iz in range(1,len(z)):
                                #i = 1D histogram bin index
                                i += 1
                                #get MR and Rsq at center of bin in 3d histogram
                                mrCenter = myTH3.GetXaxis().GetBinCenter(ix)
                                rsqCenter = myTH3.GetYaxis().GetBinCenter(iy)
                                if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=bn,
                                        sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                        referenceHist=referenceHist):
                                    #print " bin",i,"(",ix,iy,iz,") matches!"
                                    matchedAtLeastOneBin = True
                else:
                    #make a TH2Poly from each xy slice of the histogram, unroll each one and attach together
                    print "Merging bins according to custom (MR-dependent) binning"
                    i = 0
                    for iz in range(1, len(z)):
                        #one xy slice of the histogram 
                        unrollRows = unrollBins[iz-1][0]
                        unrollCols = unrollBins[iz-1][1]
                        poly = macro.makeTH2PolyFromColumns("poly"+str(iz)+name, 'poly', unrollRows, unrollCols)
                        polyBins = poly.GetBins()
                        for sigBN in range(1, poly.GetNumberOfBins()+1):
                            i += 1
                            thisSigBin = polyBins.At(sigBN-1)
                            #get MR and Rsq at center of bin in 3d histogram
                            mrCenter = (thisSigBin.GetXMax() + thisSigBin.GetXMin())/2.0
                            rsqCenter = (thisSigBin.GetYMax() + thisSigBin.GetYMin())/2.0
                            if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=bn,
                                    sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                    referenceHist=referenceHist):
                                #print " bin",i,"(",sigBN,iz,") matches!"
                                matchedAtLeastOneBin = True
                        poly.Delete()
                #don't save the histogram if there is no change from the nominal
                if not matchedAtLeastOneBin:
                    print "No matching signal bins -- discarding histogram"
                    del hists[newHistName]
        else: #TH2F case
            for bx in range(1,referenceHist.GetNbinsX()+1):
                for by in range(1,referenceHist.GetNbinsY()+1):
                    b = referenceHist.GetBin(bx,by)
                    newHistName = makeNewHistogramForUncorrelateSFs(name, centerName, systName, b, hists)
                    matchedAtLeastOneBin = False
                    #find bins in hists[name] that lie inside bin b of referenceHist
                    if unrollBins is None:
                        i = 0
                        for ix in range(1,len(x)):
                            for iy in range(1,len(y)):
                                for iz in range(1,len(z)):
                                    #i = 1D histogram bin index
                                    i+= 1
                                    #get MR and Rsq at center of bin in 3d histogram
                                    mrCenter = myTH3.GetXaxis().GetBinCenter(ix)
                                    rsqCenter = myTH3.GetYaxis().GetBinCenter(iy)
                                    if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=b,
                                            sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                            referenceHist=referenceHist):
                                        #print " bin",i,"(",ix,iy,iz,") matches!"
                                        matchedAtLeastOneBin = True
                    else:
                        #make a TH2Poly from each xy slice of the histogram, unroll each one and attach together
                        print "Merging bins according to custom (MR-dependent) binning"
                        i = 0
                        for iz in range(1, len(z)):
                            #one xy slice of the histogram 
                            unrollRows = unrollBins[iz-1][0]
                            unrollCols = unrollBins[iz-1][1]
                            poly = macro.makeTH2PolyFromColumns("poly"+str(iz)+name, 'poly', unrollRows, unrollCols)
                            polyBins = poly.GetBins()
                            for sigBN in range(1, poly.GetNumberOfBins()+1):
                                i += 1
                                thisSigBin = polyBins.At(sigBN-1)
                                #get MR and Rsq at center of bin in 3d histogram
                                mrCenter = (thisSigBin.GetXMax() + thisSigBin.GetXMin())/2.0
                                rsqCenter = (thisSigBin.GetYMax() + thisSigBin.GetYMin())/2.0
                                if setBinContentsForUncorrelateSFs(mrCenter, rsqCenter, refBN=b,
                                        sigBN=i, sysHist=hists[name], newHist=hists[newHistName],
                                        referenceHist=referenceHist):
                                    #print " bin",i,"(",sigBN,iz,") matches!"
                                    matchedAtLeastOneBin = True
                            poly.Delete()
                    #don't save the histogram if there is no change from the nominal
                    if not matchedAtLeastOneBin:
                        print "No matching signal bins -- discarding histogram"
                        del hists[newHistName]

        #remove the original histogram
        del hists[name]

def convertTree2TH1(tree, cfg, box, workspace, f, globalScaleFactor, treeName, sfs={}, sysErrOpt="", sumPdfWeights=None, sumScaleWeights=None, nevents=None, isData=False, unrollBins=None, xBR=-1, yBR=-1):
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

    #if QCD, use data-driven extrapolation
    if treeName.lower() == 'qcd':
        isData = True
        cuts = cuts.replace('abs(dPhiRazor) <','abs(dPhiRazor) >')

    #noise flags
    if isData:
        flagCuts = ' && '.join(['Flag_HBHENoiseFilter','Flag_HBHEIsoNoiseFilter','Flag_goodVertices','Flag_eeBadScFilter'])
        cuts = cuts + ' && ( ' + flagCuts + ' ) '

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
        if isData:
            theWeight = 1.0
            if treeName.lower() == 'qcd': #extrapolate from low deltaPhi region for QCD prediction
                qcdExtrapolationFactor = 1.2e+6*(tree.MR**(-2.6)) + 0.064
                theWeight *= qcdExtrapolationFactor
        else:
            theWeight = tree.weight*k*globalScaleFactor
            if xBR>-1 and yBR>-1:
                theWeight *= getBRWeight(xBR, yBR, tree.ntFromGluino, tree.nCharginoFromGluino)
                theWeight /= getBRWeight(0.25, 0.25, tree.ntFromGluino, tree.nCharginoFromGluino)
                #print getBRWeight(xBR, yBR, tree.ntFromGluino, tree.nCharginoFromGluino)/getBRWeight(0.25, 0.25, tree.ntFromGluino, tree.nCharginoFromGluino)
        filledWeight = fillRazor3D(tree, myTH3, theWeight, btagCutoff, treeName, sfs, sysErrOpt, sumPdfWeights=sumPdfWeights, sumScaleWeights=sumScaleWeights, nevents=nevents)
        numEntriesByBtag[btag_bin] += 1
        sumEntriesByBtag[btag_bin] += filledWeight

    #unroll into TH1F
    if unrollBins is None:
        nBins = (len(x)-1)*(len(y)-1)*(len(z)-1)
        maxBins = 224
        #maxBins = nBins
        myTH1 = rt.TH1F(treeName,treeName,maxBins,0,maxBins)
        myTH1.SetDirectory(0) #prevent it from going out of scope
        myTH1.Sumw2()
        i = 0
        for ix in range(1,len(x)):
            for iy in range(1,len(y)):
                for iz in range(1,len(z)):
                    i += 1
                    if sysErrOpt == "mcstatUp":
                        myTH1.SetBinContent(i,myTH3.GetBinContent(ix,iy,iz) + myTH3.GetBinError(ix,iy,iz))
                        myTH1.SetBinError(i,myTH3.GetBinError(ix,iy,iz))
                    elif sysErrOpt == "mcstatDown":
                        myTH1.SetBinContent(i,max(0.,myTH3.GetBinContent(ix,iy,iz) - myTH3.GetBinError(ix,iy,iz)))
                        myTH1.SetBinError(i,myTH3.GetBinError(ix,iy,iz))
                    elif sysErrOpt == "pdfUp": #inflate bin contents by 10%
                        myTH1.SetBinContent(i,myTH3.GetBinContent(ix,iy,iz)*1.1)
                        myTH1.SetBinError(i,myTH3.GetBinError(ix,iy,iz)*1.1)
                    elif sysErrOpt == "pdfDown": #deflate bin contents by 10%
                        myTH1.SetBinContent(i,myTH3.GetBinContent(ix,iy,iz)/1.1)
                        myTH1.SetBinError(i,myTH3.GetBinError(ix,iy,iz)/1.1)
                    else:                    
                        myTH1.SetBinContent(i,myTH3.GetBinContent(ix,iy,iz))
                        myTH1.SetBinError(i,myTH3.GetBinError(ix,iy,iz))
    else:
        #make a TH2Poly from each xy slice of the histogram, unroll each one and attach together
        print "Merging bins according to custom (MR-dependent) binning"
        layers = []
        for iz in range(1, len(z)):
            #one xy slice of the histogram (doing it myself to avoid any ROOT eccentricities)
            tempTH2 = rt.TH2D(treeName+'2d'+str(iz), 'tmp', len(x)-1,x,len(y)-1,y)
            for ix in range(1, len(x)):
                for iy in range(1, len(y)):
                    tempTH2.SetBinContent(ix,iy, myTH3.GetBinContent(ix,iy,iz))
                    tempTH2.SetBinError(ix,iy, myTH3.GetBinError(ix,iy,iz))
            #turn it into a TH2Poly with the reduced binning
            unrollRows = unrollBins[iz-1][0]
            unrollCols = unrollBins[iz-1][1]
            poly = macro.makeTH2PolyFromColumns(tempTH2.GetName()+"poly", 'poly', unrollRows, unrollCols)
            macro.fillTH2PolyFromTH2(tempTH2, poly)
            numbins = poly.GetNumberOfBins()
            unrolledSlice = rt.TH1D(tempTH2.GetName()+"Unroll", "slice", numbins, 0, numbins)
            for bn in range(1, numbins+1):
                unrolledSlice.SetBinContent(bn, poly.GetBinContent(bn))
                unrolledSlice.SetBinError(bn, poly.GetBinError(bn))
            layers.append(unrolledSlice)
            poly.Delete()
        myTH1 = macro.stitch(layers)
        myTH1.SetName(treeName)
        myTH1.SetTitle(treeName)

    print "Filename: %s"%f
    print "Sample: %s"%treeName
    print "Scale Factors     [ %s ] ="%box, k
    print "Number of Entries [ %s ] ="%(box), numEntriesByBtag
    print "Sum of Weights    [ %s ] ="%(box), sumEntriesByBtag

    return myTH1

def writeDataCard_th1(box,model,txtfileName,hists,bkgs=None):
    if bkgs is None: 
        bkgs = [bkg for bkg in FILENAMES_MC if bkg in hists]
    obsRate = hists["data_obs"].Integral()
    nBkgd = len(bkgs)
    rootFileName = txtfileName.replace('.txt','.root')
    rates = [hists[model].Integral()]
    rates.extend([hists[bkg].Integral() for bkg in bkgs])
    processes = [model]
    processes.extend(bkgs)
    lumiErrs = [1+lumiUncertainty] 
    lumiErrs.extend([1+lumiUncertainty if bkg.lower() != 'qcd' else 1.0 for bkg in bkgs]) 
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

    #20% normalization uncertainty on rare backgrounds
    for bkg in bkgs:
        if bkg in ['ttjets','wjetstolnu','dyjetstoll','zjetstonunu','qcd']: continue
        #alternate naming convention
        if bkg.lower() in ['ttjets', 'ttjets1l', 'ttjets2l', 'wjets', 'dyjets', 'zinv', 'qcd']: continue
        mcErrs[bkg] = [1.00]
        mcErrs[bkg].extend([1.00 + 0.20*(bkg==bkg1) for bkg1 in bkgs]) 
            
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
        if bkg not in mcErrs: continue
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
        if bkg in mcErrStrings: datacard+=mcErrStrings[bkg] #MC normalization uncertainties
    for name in shapeNames:
        datacard+=shapeErrStrings[name] #shape uncertainties

    #print out card
    print "\n",datacard,"\n"

    #write card
    txtfile = open(txtfileName,"w")
    txtfile.write(datacard)
    txtfile.close()
 
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--config',dest="config",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_argument('-d','--dir',dest="outDir",default="./",
                  help="Output directory to store datasets")
    parser.add_argument('-l','--lumi',dest="lumi", default=2300.,type=float,
                  help="integrated luminosity in pb^-1")
    parser.add_argument('--lumi-in',dest="lumi_in", default=1.,type=float,
                  help="integrated luminosity in pb^-1")
    parser.add_argument('-b','--box',default="MultiJet",
                  help="box name")
    parser.add_argument('-m','--model',default="T1bbbb",
                  help="signal model name")
    parser.add_argument('--mGluino',default=-1,type=float,
                  help="mass of gluino")
    parser.add_argument('--mStop',default=-1,type=float,
                  help="mass of stop")
    parser.add_argument('--mLSP',default=-1,type=float,
                  help="mass of LSP")
    parser.add_argument('--no-sys',dest="noSys",default=False,action='store_true',
                  help="no shape systematic uncertainties")
    parser.add_argument('--unblind',default=False,action='store_true',
                  help='set limits on data (default is to use MC cocktail)')

    args = parser.parse_args()
    
    cfg = Config.Config(args.config)

    box =  args.box
    boxList = box.split('_')
    lumi = args.lumi
    lumi_in = args.lumi_in

    #get signal mass point
    if args.mGluino>-1:
        massPoint = '%i_%i'%(args.mGluino,args.mLSP)
    elif args.mStop>-1:
        massPoint = '%i_%i'%(args.mStop,args.mLSP)
    modelString = 'SMS-'+args.model+'_'+massPoint

    #get theory cross section
    if args.mGluino > -1:
        thyXsec, thyXsecErr = getTheoryCrossSectionAndError(mGluino=args.mGluino)
    elif args.mStop > -1:
        thyXsec, thyXsecErr = getTheoryCrossSectionAndError(mStop=args.mStop)

    for curBox in boxList:
        #get data/MC scale factors from files
        processNames = ["ttjets", "wjetstolnu", "dyjetstoll", "zjetstonunu"]
        scaleFactorNames = {"ttjets":"TTJets","wjetstolnu":"WJets","dyjetstoll":"DYJets","zjetstonunu":"WJetsInv"}
        sfHists = loadScaleFactorHists(sfFilename="data/ScaleFactors/RazorScaleFactors_MultiJet.root", processNames=processNames, scaleFactorNames=scaleFactorNames, debugLevel=1)

        #TODO: commit cross check scale factors and access them from central location
        #cross check scale factors
        ttjetsDileptonSFHists = loadScaleFactorHists(sfFilename="RazorTTJetsDileptonCrossCheck.root", processNames=["ttjetsdilepton"], scaleFactorNames={ "ttjetsdilepton":"TTJets" }, debugLevel=1)
        vetoLeptonSFHists = loadScaleFactorHists(sfFilename="RazorVetoLeptonCrossCheck.root", processNames=["vetolepton"], scaleFactorNames={ "vetolepton":"VetoLepton" }, debugLevel=1)
        vetoTauSFHists = loadScaleFactorHists(sfFilename="RazorVetoTauCrossCheck.root", processNames=["vetotau"], scaleFactorNames={ "vetotau":"VetoTau" }, debugLevel=1)
        sfHists.update(ttjetsDileptonSFHists) #combine scale factor dictionaries
        sfHists.update(vetoLeptonSFHists) 
        sfHists.update(vetoTauSFHists) 

        #list of shape systematics to apply to signal and background MC.
        #if a list of physics processes is given, the uncertainty will be applied to each process in the list, assumed uncorrelated from process to process.
        #if an empty list is given, the uncertainty will be applied (correlated) to all processes, including signal.
        shapes = {
                'tightmuoneff':[],
                'tighteleeff':[],
                'vetomuoneff':[],
                'vetoeleeff':[],
                'jes':[],
                'muontrig':[],
                'eletrig':[],
                'btag':[],
                'facscale':[],
                'renscale':[],
                'facrenscale':[],
                'ees':[],
                'mes':[],
                'pileup':[],
                'isr':[],
                'sfstat':['ttjets','wjetstolnu','zjetstonunu'],
                'sfsys':['ttjets','wjetstolnu','zjetstonunu'],
                'ttcrosscheck':['ttjets'],
                #'vetolepcrosscheck':[],
                #'vetotaucrosscheck':[],
                'mcstat%s'%curBox.lower():([f for f in FILENAMES_MC if f.lower() != 'qcd']+['signal'])
        }
        #list of shapes that apply to signal only
        signalShapes = {
                'tightmuonfastsim':[],
                'tightelefastsim':[],
                'vetomuonfastsim':[],
                'vetoelefastsim':[],
                'btagfastsim':[],
        }
        #specify which systematics should be treated as uncorrelated bin-by-bin
        uncorrShapes = [
                "mcstat%s"%curBox.lower(),
                ] 

        #scale factor systematics are correlated according to scale factor binning
        uncorrSFShapes = [
                'sfstat',
                'ttcrosscheck',
                #'vetolepcrosscheck',
                #'vetotaucrosscheck',
                ] 

        if args.noSys:
            shapes = {}
            signalShapes = {}
            uncorrShapes = []
            uncorrSFShapes = []

        #create workspace
        w = rt.RooWorkspace("w"+curBox)
        variables = initializeWorkspace(w,cfg,curBox)    
        
        #list of histograms
        ds = []
            
        btagMin =  w.var('nBtag').getMin()
        btagMax =  w.var('nBtag').getMax()
        z = array('d', cfg.getBinning(curBox)[2]) # nBtag binning

        #make MC background histograms
        for treeName, f in FILENAMES_MC.iteritems():
            rootFile = rt.TFile.Open(f) #open file
            assert rootFile
            tree = rootFile.Get('RazorInclusive') #get tree
            #get histograms for sum of pdf and scale weights
            if treeName.lower() != 'qcd':
                if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes or 'n0pdf' in shapes:
                    nevents = rootFile.Get('NEvents')
                    assert nevents
                else:
                    nevents = None
                if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes:
                    sumScaleWeights = rootFile.Get('SumScaleWeights')
                    assert sumScaleWeights
                else:
                    sumScaleWeights = None
            else:
                nevents = None
                sumScaleWeights = None

            #add histogram to output file
            print("Building histogram for "+treeName)
            ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=treeName, sfs=sfHists))
            #get up/down histograms for shape systematics
            for shape in shapes:
                for updown in ["Up", "Down"]:
                    if shapes[shape] == [] and treeName.lower() != 'qcd':
                        print("Building histogram for "+treeName+"_"+shape+updown)
                        ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=treeName+"_"+shape+updown, sfs=sfHists, sysErrOpt=shape+updown, sumScaleWeights=sumScaleWeights, nevents=nevents))
                    elif treeName.lower() in [s.lower() for s in shapes[shape]]:
                        print("Building histogram for "+treeName+"_"+shape+(treeName.replace('_',''))+updown)
                        ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=lumi/lumi_in, treeName=treeName+"_"+shape+(treeName.replace('_',''))+updown, sfs=sfHists, sysErrOpt=shape+updown, sumScaleWeights=sumScaleWeights, nevents=nevents))
            rootFile.Close()
        #signal process
        f = DIR_SIGNAL+'/'+modelString+'.root'
        rootFile = rt.TFile.Open(f)
        assert rootFile
        tree = rootFile.Get('RazorInclusive') #get tree
        nEvents = rootFile.Get('NEvents').Integral() #get number of events processed for this mass point
        #get histograms for sum of pdf and scale weights
        if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes or 'n0pdf' in shapes:
            nevents = rootFile.Get('NEvents')
            assert nevents
        else:
            nevents = None
        if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes:
            sumScaleWeights = rootFile.Get('SumScaleWeights')
            assert sumScaleWeights
        else:
            sumScaleWeights = None
        globalScaleFactor = thyXsec*lumi/lumi_in/nEvents 
        #add histogram to output file
        print("Building histogram for "+modelString)
        ds.append(convertTree2TH1(tree, cfg, curBox, w, f , globalScaleFactor=globalScaleFactor, treeName=modelString, sfs=sfHists))
        #systematics
        allShapes = shapes.copy()
        allShapes.update(signalShapes)
        for shape in allShapes:
            for updown in ["Up", "Down"]:
                if allShapes[shape] == []:
                    print("Building histogram for "+modelString+"_"+shape+updown)
                    ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=globalScaleFactor, treeName=modelString+"_"+shape+updown, sfs=sfHists, sysErrOpt=shape+updown, sumScaleWeights=sumScaleWeights, nevents=nevents))
                elif "signal" in [s.lower() for s in allShapes[shape]]:
                    print("Building histogram for "+modelString+"_"+shape+(modelString.replace('_',''))+updown)
                    ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=globalScaleFactor, treeName=modelString+"_"+shape+"signal"+updown, sfs=sfHists, sysErrOpt=shape+updown, sumScaleWeights=sumScaleWeights, nevents=nevents))
        rootFile.Close()
        #data
        if args.unblind:
            f = DATA_NAMES[curBox]
            rootFile = rt.TFile.Open(f)
            assert rootFile
            tree = rootFile.Get('RazorInclusive') #get tree
            print "Building histogram for data"
            data = convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=1.0, treeName='data_obs', isData=True)
            rootFile.Close()
        else: #use sum of MC histograms as proxy for the data
            data = ds[0].Clone('data_obs')
            data.Reset()
            for h in ds:
                if h.GetName() in FILENAMES_MC: 
                    data = data + h

        #convert dataset list to dict
        dsDict = {}
        for d in ds: dsDict[d.GetName()] = d

        #perform uncorrelation procedure
        for shape in uncorrShapes:
            uncorrelate(dsDict, shape)
        for shape in uncorrSFShapes:
            uncorrelateSFs(dsDict, shape, sfHists, cfg, curBox)

        #output file name
        if btagMax>btagMin+1:
            outFileName = 'RazorInclusive_Histograms_lumi-%.1f_%i-%ibtag_%s.root'%(lumi/1000.,btagMin,btagMax-1,curBox)
        else:
            outFileName = 'RazorInclusive_Histograms_lumi-%.1f_%ibtag_%s.root'%(lumi/1000.,btagMin,curBox)

        #output file
        os.system('mkdir -p '+args.outDir)
        print "Output File: %s"%(args.outDir+"/"+outFileName)
        outFile = rt.TFile.Open(args.outDir+"/"+outFileName,'recreate')
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
        writeDataCard_th1(curBox,modelString,(args.outDir+"/"+outFileName).replace('.root','.txt'),dsDict)

    #run combine
    if len(boxList) == 1:
        #get card name
        if btagMax>btagMin+1:
            cardName = '%s/RazorInclusive_Histograms_lumi-%.1f_%i-%ibtag_%s.txt'%(args.outDir,lumi/1000.,btagMin,btagMax-1,boxList[0])
        else:
            cardName = '%s/RazorInclusive_Histograms_lumi-%.1f_%ibtag_%s.txt'%(args.outDir,lumi/1000.,btagMin,boxList[0])

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
        exec_me('cd '+args.outDir+'; combineCards.py '+(' '.join(cardNames))+' > RazorInclusive_Histograms_'+('_'.join(boxList))+'.txt; cd ..', False)
        #call combine
        exec_me('combine -M Asymptotic '+args.outDir+'/RazorInclusive_Histograms_'+('_'.join(boxList))+'.txt', False)


