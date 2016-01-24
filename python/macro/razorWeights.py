##Weight and trigger utilities for inclusive razor analysis

import ROOT as rt
import sys

#####################################
### WEIGHT AND TRIGGER INFO
#####################################

#QCD systematic error
QCDNORMERRFRACTION = 0.87

def getQCDExtrapolationFactor(MR):
    """Get QCD extrapolation factor as a function of MR"""
    return 3.1e+7*(MR**(-3.1)) + 0.062

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
            print "Error in pileupWeight: tree does not have a branch for",puBranch
            sys.exit()
        pileupWeight = wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(getattr(event, puBranch)))
        if debugLevel > 1: print puBranch," weight:",pileupWeight,"(",puBranch,"=",getattr(event, puBranch),")"
        return pileupWeight
    else:
        print "Error in pileupWeight: pileup reweighting histogram not found!"
        sys.exit()

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
        print "Error option:",errorOpt
        print "Weight from ntuple:",event.weight
        print "Scale by:",scale

    if len(lweightOpts) > 0:

        #QCD extrapolation weight
        if 'datadrivenqcd' in lweightOpts: 
            qcdExtrapolationFactor = getQCDExtrapolationFactor(event.MR)
            if debugLevel > 1:
                print "QCD extrapolation factor:",qcdExtrapolationFactor
            eventWeight *= qcdExtrapolationFactor

        #pileup reweighting
        if str.lower("doNPVWeights") in lweightOpts:
            eventWeight *= pileupWeight(event, wHists, debugLevel=debugLevel)
        elif str.lower("doNVtxWeights") in lweightOpts:
            eventWeight *= pileupWeight(event, wHists, puBranch="nVtx", debugLevel=debugLevel)

        #lepton scale factors
        if str.lower("doLep1Weights") in lweightOpts:
            doLep2 = (str.lower("doLep2Weights") in lweightOpts)
            eventWeight *= leptonWeight(event, wHists, doLep2, debugLevel=debugLevel)

        #trigger scale factors
        if str.lower("do1LepTrigWeights") in lweightOpts:
            doLep2Trig = (str.lower("doLep2TrigWeights") in lweightOpts)
            eventWeight *= leptonTriggerWeight(event, wHists, doLep2Trig, debugLevel=debugLevel)

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

    if 'datadrivenqcd' not in lweightOpts:
        eventWeight = scale
        if debugLevel > 1: print("Applying a weight of "+str(eventWeight))
        return eventWeight
    else: #data-driven QCD estimate
        qcdExtrapolationFactor = getQCDExtrapolationFactor(event.MR)
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

def loadScaleFactorHists(sfFilename="RazorScaleFactors.root", processNames=[], scaleFactorNames={}, debugLevel=0):
    """Returns a dict with available scale factor histograms"""
    sfFile = rt.TFile.Open(sfFilename)
    assert sfFile
    sfHists = {}
    for pname in processNames:
        histname = pname
        if pname in scaleFactorNames:
            histname = scaleFactorNames[pname]
        if debugLevel > 0: print "Looking for scale factor histogram",histname,"for",pname,"...",        

        tmp = sfFile.Get(histname+"ScaleFactors")
        if tmp: 
            sfHists[pname] = tmp
            sfHists[pname].SetDirectory(0)
            if debugLevel > 0: print "Found!"
        else:
            if debugLevel > 0: print ""

        #get up/down histograms
        tmp = sfFile.Get(histname+"ScaleFactorsUp")
        if tmp: 
            sfHists[pname+"Up"] = tmp
            sfHists[pname+"Up"].SetDirectory(0)
            if debugLevel > 0: print "Up histogram found!"
        tmp = sfFile.Get(histname+"ScaleFactorsDown")
        if tmp: 
            sfHists[pname+"Down"] = tmp
            sfHists[pname+"Down"].SetDirectory(0)
            if debugLevel > 0: print "Down histogram found!"

    sfFile.Close()
    return sfHists

def getSFHistNameForErrorOpt(errorOpt, name):
    """Returns the key in the scale factor histogram dictionary for the given error option"""
    if errorOpt in ['sfsysttjetsUp','sfsyswjetsUp','sfsyszinvUp']:
        return name+'Up'
    elif errorOpt in ['sfsysttjetsDown','sfsyswjetsDown','sfsyszinvDown']:
        return name+'Down'
    else:
        return name

def getAuxSFsForErrorOpt(auxSFs={}, errorOpt=""):
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
    elif 'sfsysvetolep' in errorOpt.lower():
        if 'VetoLepton' in auxSFs: del auxSFs['VetoLepton']
        if 'Up' in errorOpt:
            histNames.append("VetoLeptonUp")
            varNames.append("leadingGenLeptonPt")
            cuts.append("abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13")
        elif 'Down' in errorOpt:
            histNames.append("VetoLeptonDown")
            varNames.append("leadingGenLeptonPt")
            cuts.append("abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13")
    #Veto tau scale factors up/down
    elif 'sfsysvetotau' in errorOpt.lower():
        if 'VetoTau' in auxSFs: del auxSFs['VetoTau']
        if 'Up' in errorOpt:
            histNames.append("VetoTauUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoTauDown")
        varNames.append("leadingGenLeptonPt")
        cuts.append("abs(leadingGenLeptonType) == 15")
    #MT efficiency up/down
    elif 'mteff' in errorOpt.lower():
        if 'VetoLepton' in auxSFs: del auxSFs['VetoLepton']
        if 'VetoTau' in auxSFs: del auxSFs['VetoTau']
        if 'Up' in errorOpt:
            histNames.append("VetoLeptonMTUp")
            histNames.append("VetoTauMTUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoLeptonMTDown")
            histNames.append("VetoTauMTDown")
        varNames.append("leadingGenLeptonPt")
        cuts.append("abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13")
        varNames.append("leadingGenLeptonPt")
        cuts.append("abs(leadingGenLeptonType) == 15")
    #DPhi efficiency up/down
    elif 'dphieff' in errorOpt.lower():
        if 'VetoLepton' in auxSFs: del auxSFs['VetoLepton']
        if 'VetoTau' in auxSFs: del auxSFs['VetoTau']
        if 'Up' in errorOpt:
            histNames.append("VetoLeptonDPhiUp")
            histNames.append("VetoTauDPhiUp")
        elif 'Down' in errorOpt:
            histNames.append("VetoLeptonDPhiDown")
            histNames.append("VetoTauDPhiDown")
        varNames.append("leadingGenLeptonPt")
        cuts.append("abs(leadingGenLeptonType) == 11 || abs(leadingGenLeptonType) == 13")
        varNames.append("leadingGenLeptonPt")
        cuts.append("abs(leadingGenLeptonType) == 15")
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
    auxSFs.update(sfsNeeded)

def invertHistogram(hist):
    """Replaces contents of each hist bin with 1/(contents).  Updates bin errors accordingly.
       For bins with no contents, does nothing."""
    ret = hist.Clone(hist.GetName()+"Inverted")
    for b in range(hist.GetSize()+1):
        if hist.GetBinContent(b) != 0:
            ret.SetBinContent( b, 1.0/hist.GetBinContent(b) )
            ret.SetBinError( b, hist.GetBinError(b) / (hist.GetBinContent(b))**2 )
    return ret

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
        'sfsysvetolep':True,
        'sfsysvetotau':True,
        'mteff':True,
        'dphieff':True,
        'singletopnorm':True,
        'othernorm':True,
        'qcdnorm':True,
        'sfsysttjets':True,
        'sfsyswjets':True,
        'sfsyszinv':True,
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
