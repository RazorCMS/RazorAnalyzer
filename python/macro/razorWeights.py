##Weight and trigger utilities for inclusive razor analysis

import ROOT as rt
import sys

#####################################
### WEIGHT AND TRIGGER INFO
#####################################

#QCD systematic error
QCDNORMERRFRACTION_DIJET = 0.85
QCDNORMERRFRACTION_MULTIJET = 0.80
#QCDNORMERRFRACTION = 0.87 #used in 2015

def getQCDExtrapolationFactor(MR,region='multijet'):
    """Get QCD extrapolation factor as a function of MR"""
    #return 3.1e+7*(MR**(-3.1)) + 0.062 #power law + constant (MultiJet 2015)
    if region.lower() == 'dijet':
        return 0.103
    else:
        return 0.242

WEIGHTDIR_DEFAULT = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors"
LEPTONWEIGHTDIR_DEFAULT = "LeptonEfficiencies/20151013_PR_2015D_Golden_1264"
weightfilenames_DEFAULT = {
        "muon": WEIGHTDIR_DEFAULT+"/"+LEPTONWEIGHTDIR_DEFAULT+"/efficiency_results_TightMuonSelectionEffDenominatorReco_2015D_Golden.root",
        "ele": WEIGHTDIR_DEFAULT+"/"+LEPTONWEIGHTDIR_DEFAULT+"/efficiency_results_TightElectronSelectionEffDenominatorReco_2015D_Golden.root",
        "muontrig": WEIGHTDIR_DEFAULT+"/"+LEPTONWEIGHTDIR_DEFAULT+"/efficiency_results_MuTriggerIsoMu27ORMu50EffDenominatorTight_2015D_Golden.root",
        "eletrig": WEIGHTDIR_DEFAULT+"/"+LEPTONWEIGHTDIR_DEFAULT+"/efficiency_results_EleTriggerEleCombinedEffDenominatorTight_2015D_Golden.root",
        "pileup": WEIGHTDIR_DEFAULT+"/PileupWeights/NVtxReweight_ZToMuMu_2015D_1264ipb.root",
        }
weighthistnames_DEFAULT = {
        "muon": "ScaleFactor_TightMuonSelectionEffDenominatorReco",
        "ele": "ScaleFactor_TightElectronSelectionEffDenominatorReco",
        "muontrig": "ScaleFactor_MuTriggerIsoMu27ORMu50EffDenominatorTight",
        "eletrig": "ScaleFactor_EleTriggerEleCombinedEffDenominatorTight",
        "pileup": "NVtxReweight",
        }

def passTrigger(event, triggerNumList):
    """Checks if the event passed any trigger in the list"""
    passes = False
    for i in triggerNumList:
        if event.HLTDecision[i]:
            passes = True
            break
    return passes

def pileupWeight(event, wHists, puBranch="NPV", debugLevel=0):
    if "pileup" in wHists:
        if not hasattr(event, puBranch):
            print "Error in pileupWeight: tree does not have a branch for ",puBranch
            sys.exit()
        pileupWeight = wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(getattr(event, puBranch)))
        if debugLevel > 1: print puBranch," weight:",pileupWeight,"(",puBranch,"=",getattr(event, puBranch),")"
        return pileupWeight
    else:
        print "Error in pileupWeight: pileup reweighting histogram not found!"
        sys.exit()

def reapplyPileupWeight(event, wHists, weightBranch="pileupWeight", debugLevel=0):
    """Get (pileupweight)/(oldpileupweight)"""
    if "pileup" in wHists:
        #use NPU_0 or NPU, whichever exists in the tree
        if not hasattr(event, "NPU_0"):
            if not hasattr(event, "NPU"):
                sys.exit( "Error in pileupWeight: tree does not have a branch for NPU_0 or NPU" )
            else:
                puBranch = "NPU"
        else:
            puBranch = "NPU_0"
        if not hasattr(event, weightBranch):
            sys.exit( "Error in pileupWeight: tree does not have a branch for "+weightBranch )
        pileupWeight = wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(
            getattr(event, puBranch))) / getattr(event,weightBranch)
        if debugLevel > 1: print puBranch," weight:",pileupWeight,"(",puBranch,"=",getattr(event, puBranch),")"
        return pileupWeight
    else:
        print "Error in pileupWeight: pileup reweighting histogram not found!"
        sys.exit()

def reapplyLepEffWeight(event, wHists, eleBranch="eleEffWeight", muBranch="muonEffWeight", debugLevel=0):
    """Get (new lepton weight)/(old lepton weight)"""
    if abs(event.lep1Type) == 13:
        if "muoneff" in wHists:
            if not hasattr(event, muBranch):
                sys.exit( "Error in reapplyLepEffWeight: tree does not have a branch for "+muBranch )
            muWeight = wHists["muoneff"].GetBinContent(
                    wHists["muoneff"].GetXaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ),
                    wHists["muoneff"].GetYaxis().FindFixBin( abs( event.lep1.Eta() ) ) ) / getattr( event, muBranch )
            if debugLevel > 1: 
                print "Muon weight:",muWeight,"(",muBranch,"=",getattr(event, muBranch),")"
            return muWeight
        else:
            print "Error: muon efficiency reweighting histogram not found!"
            sys.exit()
    elif abs(event.lep1Type) == 11:
        if "eleeff" in wHists:
            if not hasattr(event, eleBranch):
                sys.exit( "Error in reapplyLepEffWeight: tree does not have a branch for "+eleBranch )
            eleWeight = wHists["eleeff"].GetBinContent(
                    wHists["eleeff"].GetXaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ),
                    wHists["eleeff"].GetYaxis().FindFixBin( abs( event.lep1.Eta() ) ) ) / getattr( event, eleBranch )
            if debugLevel > 1: 
                print "Electron weight:",eleWeight,"(",eleBranch,"=",getattr(event, eleBranch),")"
            return eleWeight
        else:
            print "Error: electron efficiency reweighting histogram not found!"
            sys.exit()

def reapplyLepTrigWeight(event, wHists, trigBranch="trigWeight1L", debugLevel=0):
    """Get (new trig weight)/(old trig weight)"""
    if not hasattr(event, trigBranch):
        sys.exit( "Error in reapplyLepTrigWeight: tree does not have a branch for "+trigBranch )
    oldWeight = getattr( event, trigBranch )
    if abs(event.lep1Type) == 13:
        if "muontrig" in wHists:
            newWeight = wHists["muontrig"].GetBinContent(
                    wHists["muontrig"].GetXaxis().FindFixBin( abs( event.lep1.Eta() ) ),
                    wHists["muontrig"].GetYaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ) ) 
        else:
            print "Error: muon trigger reweighting histogram not found!"
            sys.exit()
    elif abs(event.lep1Type) == 11:
        if "eletrig" in wHists:
            newWeight = wHists["eletrig"].GetBinContent(
                    wHists["eletrig"].GetXaxis().FindFixBin( abs( event.lep1.Eta() ) ),
                    wHists["eletrig"].GetYaxis().FindFixBin( max( min(event.lep1.Pt(), 199.9), 20.01 ) ) )
        else:
            print "Error: electron trigger reweighting histogram not found!"
            sys.exit()
    if debugLevel > 1: 
        print "Trigger weight:",newWeight,"( old weight =",oldWeight,")"
    return newWeight / oldWeight

def getNBJetsWeight(event, debugLevel=0):
    """Reweight according to number of gen-level b-jets"""
    sfTag = 0.95
    sfNoTag = (1/.68 - sfTag) / (1/.68 - 1) #using 0.68 as the average b-tag efficiency
    bWeight = 1.0
    nTags = event.NBJetsMedium
    nGenB = event.NGenBJets
    for nb in range(nGenB):
        if nTags > nb: #we successfully tagged this many jets
            bWeight *= sfTag
        else: #we didn't tag this many jets
            bWeight *= sfNoTag
    if debugLevel > 1:
        print "Found",nGenB,"gen b-jets in event and tagged",nTags,"; assigning a weight of",bWeight
    return bWeight

def leptonWeight(event, wHists, doLep2=False, debugLevel=0):
    """Weights using leading lepton pt and eta.  If doLep2 is True, weights using subleading lepton too."""
    #check for histograms
    if not ("ele" in wHists):
        print "Error in leptonWeight: electron scale factor histogram not found!"
        sys.exit()
    if not ("muon" in wHists):
        print "Error in leptonWeight: muon scale factor histogram not found!"
        sys.exit()

    leptonWeight = 1.0
    #leading muon 
    if abs(event.lep1Type) == 13 and event.lep1.Pt() > 0:
        leptonWeight *= wHists["muon"].GetBinContent(
                wHists["muon"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["muon"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["muon"].GetXaxis().GetXmax()-1.0)),
                wHists["muon"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    #leading electron
    elif abs(event.lep1Type) == 11 and event.lep1.Pt() > 0:
        leptonWeight *= wHists["ele"].GetBinContent(
                wHists["ele"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["ele"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["ele"].GetXaxis().GetXmax()-1.0)),
                wHists["ele"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    if doLep2:
        #subleading muon 
        if abs(event.lep2Type) == 13 and event.lep2.Pt() > 0:
            leptonWeight *= wHists["muon"].GetBinContent(
                    wHists["muon"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["muon"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["muon"].GetXaxis().GetXmax()-1.0)),
                    wHists["muon"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
        #subleading electron
        elif abs(event.lep2Type) == 11 and event.lep2.Pt() > 0:
            leptonWeight *= wHists["ele"].GetBinContent(
                    wHists["ele"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["ele"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["ele"].GetXaxis().GetXmax()-1.0)),
                    wHists["ele"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
    if debugLevel > 1: 
        print "lepton weight:",leptonWeight
    return leptonWeight

def leptonTriggerWeight(event, wHists, doLep2=False, debugLevel=0):
    if not ("eletrig" in wHists):
        print "Error in leptonTriggerWeight: electron trigger scale factor histogram not found!"
        sys.exit()
    if not ("muontrig" in wHists):
        print "Error in leptonTriggerWeight: muon trigger scale factor histogram not found!"
        sys.exit()
    trigWeight = 1.0
    #leading muon 
    if abs(event.lep1Type) == 13 and event.lep1.Pt() > 0:
        trigWeight *= wHists["muontrig"].GetBinContent(
                wHists["muontrig"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["muontrig"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["muontrig"].GetXaxis().GetXmax()-1.0)),
                wHists["muontrig"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    #leading electron
    elif abs(event.lep1Type) == 11 and event.lep1.Pt() > 0:
        trigWeight *= wHists["eletrig"].GetBinContent(
                wHists["eletrig"].GetXaxis().FindFixBin(max(min(event.lep1.Pt(), 199.9),15.01)),
                #wHists["eletrig"].GetXaxis().FindFixBin(min(event.lep1.Pt(), wHists["eletrig"].GetXaxis().GetXmax()-1.0)),
                wHists["eletrig"].GetYaxis().FindFixBin(abs(event.lep1.Eta())))
    if doLep2: 
        #subleading muon 
        if abs(event.lep2Type) == 13 and event.lep2.Pt() > 0:
            trigWeight *= wHists["muontrig"].GetBinContent(
                    wHists["muontrig"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["muontrig"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["muontrig"].GetXaxis().GetXmax()-1.0)),
                    wHists["muontrig"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
        #subleading electron
        elif abs(event.lep2Type) == 11 and event.lep2.Pt() > 0:
            trigWeight *= wHists["eletrig"].GetBinContent(
                    wHists["eletrig"].GetXaxis().FindFixBin(max(min(event.lep2.Pt(), 199.9),15.01)),
                    #wHists["eletrig"].GetXaxis().FindFixBin(min(event.lep2.Pt(), wHists["eletrig"].GetXaxis().GetXmax()-1.0)),
                    wHists["eletrig"].GetYaxis().FindFixBin(abs(event.lep2.Eta())))
    if debugLevel > 1: 
        if not doLep2:
            print "1-lepton trigger weight:",trigWeight
        else: 
            print "2-lepton trigger weight:",trigWeight
    return trigWeight

def weight_mc(event, wHists, scale=1.0, weightOpts=[], errorOpt=None, debugLevel=0):
    """Apply pileup weights and other known MC correction factors"""
    lweightOpts = map(str.lower, weightOpts)

    #nominal weight
    eventWeight = event.weight*scale
    if debugLevel > 1: 
        print "Weight options:",weightOpts
        print "Error option:",errorOpt
        print "Weight from ntuple:",event.weight
        print "Scale by:",scale

    if len(lweightOpts) > 0:

        #QCD extrapolation weight
        if 'datadrivenqcddijet' in lweightOpts: 
            qcdExtrapolationFactor = getQCDExtrapolationFactor(event.MR,region='dijet')
            if debugLevel > 1:
                print "QCD extrapolation factor:",qcdExtrapolationFactor
            eventWeight *= qcdExtrapolationFactor
        elif 'datadrivenqcdmultijet' in lweightOpts: 
            qcdExtrapolationFactor = getQCDExtrapolationFactor(event.MR,region='multijet')
            if debugLevel > 1:
                print "QCD extrapolation factor:",qcdExtrapolationFactor
            eventWeight *= qcdExtrapolationFactor
        elif 'qcdphoton' in lweightOpts:
            qcdWeight = 0.05
            if debugLevel > 1:
                print "QCD weight:",qcdWeight
            eventWeight *= qcdWeight 

        #reweighting in number of b-jets
        if 'nbjets' in lweightOpts:
            eventWeight *= getNBJetsWeight(event, debugLevel=debugLevel)

        #top pt reweighting
        if 'toppt' in lweightOpts and hasattr(event, 'topPtWeight') and event.topPtWeight > 0:
            if debugLevel > 1:
                print "Top pt weight:",event.topPtWeight
            eventWeight *=  event.topPtWeight

        #pileup reweighting
        if str.lower("reapplyNPUWeights") in lweightOpts:
            eventWeight *= reapplyPileupWeight(event, wHists, debugLevel=debugLevel)
        #elif str.lower("doNPVWeights") in lweightOpts:
        #    eventWeight *= pileupWeight(event, wHists, debugLevel=debugLevel)
        #elif str.lower("doNVtxWeights") in lweightOpts:
        #    eventWeight *= pileupWeight(event, wHists, puBranch="nVtx", debugLevel=debugLevel)

        #lepton scale factors
        if str.lower("reapplyLepWeights") in lweightOpts:
            eventWeight *= reapplyLepEffWeight(event, wHists, debugLevel=debugLevel)
        if str.lower("reapplyTrigWeights") in lweightOpts:
            eventWeight *= reapplyLepTrigWeight(event, wHists, debugLevel=debugLevel)
        elif str.lower("removeTrigWeights") in lweightOpts:
            eventWeight /= event.trigWeight1L
        if str.lower("removePileupWeights") in lweightOpts:
            eventWeight /= event.pileupWeight
        #if str.lower("doLep1Weights") in lweightOpts:
        #    doLep2 = (str.lower("doLep2Weights") in lweightOpts)
        #    eventWeight *= leptonWeight(event, wHists, doLep2, debugLevel=debugLevel)

        #trigger scale factors
        #if str.lower("do1LepTrigWeights") in lweightOpts:
        #    doLep2Trig = (str.lower("doLep2TrigWeights") in lweightOpts)
        #    eventWeight *= leptonTriggerWeight(event, wHists, doLep2Trig, debugLevel=debugLevel)

        if 'photonkfactor' in lweightOpts:
            kfactor = 1.44
            if debugLevel > 1:
                print "Photon k-factor",kfactor
            eventWeight *= kfactor

    #up/down corrections for systematics
    normErrFraction=0.2
    if errorOpt is not None:
        if errorOpt == "tightmuoneffUp":
            eventWeight *= event.sf_muonEffUp
            if debugLevel > 1: print "muonEffUp scale factor:",event.sf_muonEffUp
        elif errorOpt == "tightmuoneffDown":
            eventWeight *= event.sf_muonEffDown
            if debugLevel > 1: print "muonEffDown scale factor:",event.sf_muonEffDown
        elif errorOpt == "tighteleeffUp":
            eventWeight *= event.sf_eleEffUp
            if debugLevel > 1: print "eleEffUp scale factor:",event.sf_eleEffUp
        elif errorOpt == "tighteleeffDown":
            eventWeight *= event.sf_eleEffDown
            if debugLevel > 1: print "eleEffDown scale factor:",event.sf_eleEffDown
        elif errorOpt == "vetomuoneffUp":
            eventWeight *= event.sf_vetoMuonEffUp
            if debugLevel > 1: print "vetoMuonEffUp scale factor:",event.sf_vetoMuonEffUp
        elif errorOpt == "vetomuoneffDown":
            eventWeight *= event.sf_vetoMuonEffDown
            if debugLevel > 1: print "vetoMuonEffDown scale factor:",event.sf_vetoMuonEffDown
        elif errorOpt == "vetoeleeffUp":
            eventWeight *= event.sf_vetoEleEffUp
            if debugLevel > 1: print "vetoEleEffUp scale factor:",event.sf_vetoEleEffUp
        elif errorOpt == "vetoeleeffDown":
            eventWeight *= event.sf_vetoEleEffDown
            if debugLevel > 1: print "vetoEleEffDown scale factor:",event.sf_vetoEleEffDown
        elif errorOpt == "muontrigUp":
            eventWeight *= event.sf_muonTrigUp
            if debugLevel > 1: print "muontrigUp scale factor:",event.sf_muonTrigUp
        elif errorOpt == "muontrigDown":
            eventWeight *= event.sf_muonTrigDown
            if debugLevel > 1: print "muontrigDown scale factor:",event.sf_muonTrigDown
        elif errorOpt == "eletrigUp":
            eventWeight *= event.sf_eleTrigUp
            if debugLevel > 1: print "eletrigUp scale factor:",event.sf_eleTrigUp
        elif errorOpt == "eletrigDown":
            eventWeight *= event.sf_eleTrigDown
            if debugLevel > 1: print "eletrigDown scale factor:",event.sf_eleTrigDown
        elif errorOpt == "btagUp":
            eventWeight *= event.sf_btagUp
            if debugLevel > 1: print "btagUp scale factor:",event.sf_btagUp
        elif errorOpt == "btagDown":
            eventWeight *= event.sf_btagDown
            if debugLevel > 1: print "btagDown scale factor:",event.sf_btagDown
        elif errorOpt == "bmistagUp":
            eventWeight *= event.sf_bmistagUp
            if debugLevel > 1: print "bmistagUp scale factor:",event.sf_bmistagUp
        elif errorOpt == "bmistagDown":
            eventWeight *= event.sf_bmistagDown
            if debugLevel > 1: print "bmistagDown scale factor:",event.sf_bmistagDown
        elif errorOpt == "pileupUp":
            eventWeight *= event.pileupWeightUp
            if debugLevel > 1: print "pileupWeightUp scale factor:",event.pileupWeightUp
        elif errorOpt == "pileupDown":
            eventWeight *= event.pileupWeightDown
            if debugLevel > 1: print "pileupWeightDown scale factor:",event.pileupWeightDown
        elif errorOpt == "isrUp":
            eventWeight *= event.ISRSystWeightUp
            if debugLevel > 1: print "ISRSystWeightUp scale factor:",event.ISRSystWeightUp
        elif errorOpt == "isrDown":
            eventWeight *= event.ISRSystWeightDown
            if debugLevel > 1: print "ISRSystWeightDown scale factor:",event.ISRSystWeightDown
        elif errorOpt == "facscaleUp":
            eventWeight *= event.sf_facScaleUp
            if debugLevel > 1: print "facScaleUp scale factor:",event.sf_facScaleUp
        elif errorOpt == "facscaleDown":
            eventWeight *= event.sf_facScaleDown
            if debugLevel > 1: print "facScaleDown scale factor:",event.sf_facScaleDown
        elif errorOpt == "renscaleUp":
            eventWeight *= event.sf_renScaleUp
            if debugLevel > 1: print "renScaleUp scale factor:",event.sf_renScaleUp
        elif errorOpt == "renscaleDown":
            eventWeight *= event.sf_renScaleDown
            if debugLevel > 1: print "renScaleDown scale factor:",event.sf_renScaleDown
        elif errorOpt == "facrenscaleUp":
            eventWeight *= event.sf_facRenScaleUp
            if debugLevel > 1: print "facRenScaleUp scale factor:",event.sf_facRenScaleUp
        elif errorOpt == "facrenscaleDown":
            eventWeight *= event.sf_facRenScaleDown
            if debugLevel > 1: print "facRenScaleDown scale factor:",event.sf_facRenScaleDown
        elif 'normUp' in errorOpt:
            eventWeight *= (1+normErrFraction)
            if debugLevel > 1: print errorOpt,"scale factor:",1+normErrFraction
        elif 'normDown' in errorOpt:
            eventWeight /= (1+normErrFraction)
            if debugLevel > 1: print errorOpt,"scale factor:",1/(1+normErrFraction)

    if debugLevel > 1: 
        print "event weight:",eventWeight
    return eventWeight

def weight_data(event, wHists, scale=1.0, weightOpts=[], errorOpt=None, debugLevel=0):
    lweightOpts = map(str.lower, weightOpts)

    if len(lweightOpts) == 0:
        eventWeight = scale
        if debugLevel > 1: print("Applying a weight of "+str(eventWeight))
        return eventWeight
    else: #data-driven QCD estimate
        if 'datadrivenqcddijet' in lweightOpts:
            qcdExtrapolationFactor = getQCDExtrapolationFactor(event.MR,region='dijet')
        elif 'datadrivenqcdmultijet' in lweightOpts:
            qcdExtrapolationFactor = getQCDExtrapolationFactor(event.MR,region='multijet')
        elif 'qcdphoton' in lweightOpts:
            qcdExtrapolationFactor = 0.05
        else:
            qcdExtrapolationFactor = 1.0
            print "Warning: data weight options",lweightOpts,"may not make sense; see macro.razorWeights.weight_data"
        eventWeight = qcdExtrapolationFactor*scale
        if debugLevel > 1:
            print "QCD extrapolation factor:",qcdExtrapolationFactor
            print "Scale by:",scale
            print "event weight:",eventWeight
        return eventWeight

def getMTRelUncertainty(MR, bkg, box):
    muonBoxes = ["MuMultiJet", "MuSixJet", "MuFourJet", "MuJet"]
    eleBoxes = ["EleMultiJet", "EleSixJet", "EleFourJet", "EleJet"]
    bkg = bkg.lower()
    unc = 0.0
    if bkg == "ttjets" or bkg == "ttbar":
        if box in muonBoxes:
            if MR < 400:
                unc = 0.0614822 #inclusive value -- no uncertainty is available for MR < 400
            elif MR < 600:
                unc = 0.0843497
            elif MR < 800:
                unc = 0.0977245
            elif MR < 1000:
                unc = 0.124865
            else:
                unc = 0.0950384
        if box in eleBoxes:
            if MR < 400:
                unc = 0.0337985
            elif MR < 600:
                unc = 0.0254958
            elif MR < 800:
                unc = 0.0338464
            elif MR < 1000:
                unc = 0.0376923
            else:
                unc = 0.0367005
    elif bkg == "wjets":
        if box in muonBoxes:
            if MR < 400:
                unc = 0.172134
            elif MR < 600:
                unc = 0.210269
            elif MR < 800:
                unc = 0.159052
            elif MR < 1000:
                unc = 0.242155
            else:
                unc = 0.263298
        if box in eleBoxes:
            if MR < 400:
                unc = 0.147676
            elif MR < 600:
                unc = 0.0656713
            elif MR < 800:
                unc = 0.0871765
            elif MR < 1000:
                unc = 0.0979534
            else:
                unc = 0.164399
    return unc

def applyMTUncertainty1D(hist, process, debugLevel=0):
    """hist is assumed to be a histogram of MR"""
    if process == "": return 

    bkg = process.split('_')[0].lower()
    box = process.split('_')[1]
    if bkg != "ttjets" and bkg != "ttbar" and bkg != "wjets": return

    #propagate MT uncertainty to each bin
    for bx in range(1,hist.GetNbinsX()+1):
        MRAtBinCenter = hist.GetXaxis().GetBinCenter(bx)
        mtRelUnc = getMTRelUncertainty(MRAtBinCenter, bkg, box)
        mtUnc = mtRelUnc*hist.GetBinContent(bx) #convert relative error to absolute
        if debugLevel > 0:
            print "MT uncertainty on bin",bx,"(",bkg,box,") is",mtUnc,("(%1.3f of %1.3f)" % (mtRelUnc,hist.GetBinContent(bx)))
        currUnc = hist.GetBinError(bx)
        hist.SetBinError(bx, (mtUnc*mtUnc + currUnc*currUnc)**(0.5))
        if debugLevel > 0:
            print "Uncertainty on this bin increases from",currUnc,"to",hist.GetBinError(bx)

def applyMTUncertainty2D(hist, process, debugLevel=0):
    """hist is assumed to be a histogram with MR on the x-axis"""
    if process == "": return 

    bkg = process.split('_')[0].lower()
    box = process.split('_')[1]
    if bkg != "ttjets" and bkg != "ttbar" and bkg != "wjets": return

    #propagate MT uncertainty to each bin
    for bx in range(1,hist.GetNbinsX()+1):
        MRAtBinCenter = hist.GetXaxis().GetBinCenter(bx)
        mtRelUnc = getMTRelUncertainty(MRAtBinCenter, bkg, box)
        for by in range(1,hist.GetNbinsY()+1):
            mtUnc = mtRelUnc*hist.GetBinContent(bx,by) #convert relative error to absolute
            if debugLevel > 0:
                print "MT uncertainty on bin",bx,",",by,"(",bkg,box,") is",mtUnc,("(%1.3f of %1.3f)" % (mtRelUnc,hist.GetBinContent(bx,by)))

            currUnc = hist.GetBinError(bx,by)
            hist.SetBinError(bx,by, (mtUnc*mtUnc + currUnc*currUnc)**(0.5))
            if debugLevel > 0:
                print "Uncertainty on this bin increases from",currUnc,"to",hist.GetBinError(bx,by)

def loadWeightHists(filenames={}, histnames={}, debugLevel=0):
    """Returns a dict with necessary reweighting histograms"""
    if "pileup" not in filenames:
        filenames["pileup"] = "data/ScaleFactors/Placeholders/DummyRun2PileupWeights.root"
        histnames["pileup"] = "PUWeight_Run2"
    if "muon" not in filenames:
        filenames["muon"] = "data/ScaleFactors/Placeholders/DummyRun2MuonWeights.root" #dummy file
        histnames["muon"] = "MuonWeight_Run2_Tight" 
    if "ele" not in filenames:
        filenames["ele"] = "data/ScaleFactors/Placeholders/DummyRun2EleWeights.root" #dummy file
        histnames["ele"] = "EleWeight_Run2_Tight" #razor tag&probe
    wHists = {}
    wFiles = {}
    for name in filenames:
        if debugLevel > 0: print "Opening",name,"weight file",filenames[name]
        wFiles[name] = rt.TFile.Open(filenames[name])
        if debugLevel > 0: print "Getting histogram",histnames[name]
        wHists[name] = wFiles[name].Get(histnames[name])
        wHists[name].SetDirectory(0)
        assert wHists[name]
    for name in wFiles:
        wFiles[name].Close()
    return wHists

def getSFHistNameForErrorOpt(errorOpt, name):
    """Returns the key in the scale factor histogram dictionary for the given error option"""
    if errorOpt in ['sfsysttjetsUp','sfsyswjetsUp','sfsyszinvUp']:
        return name+'Up'
    elif errorOpt in ['sfsysttjetsDown','sfsyswjetsDown','sfsyszinvDown']:
        return name+'Down'
    else:
        return name

vetoLeptonAuxCuts="(abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13) && leadingGenLeptonPt > 5"
vetoTauAuxCuts="abs(leadingGenLeptonType) == 15 && leadingGenLeptonPt > 20"
def getAuxSFsForErrorOpt(auxSFs={}, errorOpt="", auxSFsPerProcess=False):
    """
    Returns scale factor histogram names needed to compute the indicated shape uncertainty.
    Format of the input/output is { "HistogramName":("variableName", "cuts"), ... }
    """

    #for building output
    histNames=[]
    varNames=[]
    cuts=[]

    ###supported error options
    #TTJets dilepton control region systematic
    if 'ttcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("TTJetsDileptonUp")
        elif 'Down' in errorOpt:
            histNames.append("TTJetsDileptonDown")
        varNames.append(("MR","Rsq"))
        cuts.append("1")
    #DYJets dilepton control region systematic
    elif 'zllcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("DYJetsInvUp")
        elif 'Down' in errorOpt:
            histNames.append("DYJetsInvDown")
        varNames.append(("MR","Rsq"))
        cuts.append("1")
    #Veto lepton scale factors up/down
    elif 'vetolepptcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoLeptonPtUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoLeptonPtDown")
        varNames.append("leadingGenLeptonPt")
        cuts.append(vetoLeptonAuxCuts)
    elif 'vetolepetacrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoLeptonEtaUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoLeptonEtaDown")
        varNames.append("abs(leadingGenLeptonEta)")
        cuts.append(vetoLeptonAuxCuts)
    #Veto tau scale factors up/down
    elif 'vetotauptcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoTauPtUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoTauPtDown")
        varNames.append("leadingGenLeptonPt")
        cuts.append(vetoTauAuxCuts)
    elif 'vetotauetacrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("VetoTauEtaUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoTauEtaDown")
        varNames.append("abs(leadingGenLeptonEta)")
        cuts.append(vetoTauAuxCuts)
    #b-tag bins closure test systematic
    elif 'btagcrosscheckmr' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("MRBUp")
        elif 'Down' in errorOpt:
            histNames.append("MRBDown")
        varNames.append("MR")
        cuts.append("1")
    elif 'btagcrosscheckrsq' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append("RsqBUp")
        elif 'Down' in errorOpt:
            histNames.append("RsqBDown")
        varNames.append("Rsq")
        cuts.append("1")
    #b-tag closure test systematic for ZInv
    elif 'btaginvcrosscheck' in errorOpt.lower():
        if 'Up' in errorOpt:
            histNames.append('ZInvBUp')
        elif 'Down' in errorOpt:
            histNames.append('ZInvBDown')
        varNames.append('nBTaggedJets')
        cuts.append('1')

    #return dictionary with needed information
    sfsNeeded = { histNames[i]:(varNames[i],cuts[i]) for i in range(len(histNames)) }
    if auxSFsPerProcess:
        for process in auxSFs:
            auxSFs[process].update(sfsNeeded.copy())
    else:
        auxSFs.update(sfsNeeded)

def splitShapeErrorsByType(shapeErrors):
    """Takes a list of shape uncertainties and splits it into two lists: the first is the list of uncertainties applied as per-event scale factors, and the second is the list of uncertainties that require separate processing."""
    supportedShapeUncertainties = { #True: belongs in list 1.  False: belongs in list 2
        'jes':False,
        'ees':False,
        'mes':False,
        'btag':True,
        'pileup':True,
        'bmistag':True,
        'facscale':True,
        'renscale':True,
        'facrenscale':True,
        'ttcrosscheck':True,
        'zllcrosscheck':True,
        'btagcrosscheckmr':True,
        'btagcrosscheckrsq':True,
        'btaginvcrosscheck':True,
        'vetolepptcrosscheck':True,
        'vetotauptcrosscheck':True,
        'vetolepetacrosscheck':True,
        'vetotauetacrosscheck':True,
        'singletopnorm':True,
        'othernorm':True,
        'qcdnorm':True,
        'sfsysttjets':True,
        'sfsyswjets':True,
        'sfsyszinv':True,
        'sfstatttjets':True,
        'sfstatwjets':True,
        'sfstatzinv':True,
        'vetomuoneff':True,
        'vetoeleeff':True,
        'tightmuoneff':True,
        'tighteleeff':True,
        'muontrig':True,
        'eletrig':True,
        }
    sfUncertainties = []
    otherUncertainties = []
    for shape in shapeErrors:
        #if shape is wrapped in a tuple, unwrap it
        if not isinstance(shape, basestring):
            thisShape = shape[0]
        else:
            thisShape = shape
        if thisShape in supportedShapeUncertainties:
            if supportedShapeUncertainties[thisShape]:
                sfUncertainties.append(shape)
            else:
                otherUncertainties.append(shape)
        else:
            print "Warning in splitShapeErrorsByType: error option",thisShape,"is not supported!"

    return sfUncertainties, otherUncertainties

def getNJetsSFs(analysis,jetName='NJets40'):
    """From an Analysis object, get the needed NJets scale factors"""
    auxSFs = { process:{} for process in analysis.samples }
    for name in ['WJets','TTJets','TTJets1L','TTJets2L']:
        if name in analysis.samples:
            auxSFs[name] = {'NJets':(jetName,'1')}
    for name in ['WJetsInv','DYJetsInv','ZInv','GJetsInv']:
        if name in analysis.samples:
            auxSFs[name] = {"NJetsInv":(jetName,"1")}
    return auxSFs

