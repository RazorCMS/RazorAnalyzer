## Inclusive razor analysis definitions
import ROOT as rt

#Single lepton trigger indices
singleLeptonTriggerNumsData = [1,2,8,20,22,24,25]
singleLeptonTriggerNumsMC = [1,2,8,17,18,19,24,25]

#TTJets Single Lepton cuts
ttjetsSingleLeptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1Pt > 30 && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium > 0 && MR > 300 && Rsq > 0.15"
ttjetsSingleLeptonCutsData = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in singleLeptonTriggerNumsData])) + ") && " + ttjetsSingleLeptonCuts
ttjetsSingleLeptonCutsMC = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in singleLeptonTriggerNumsMC])) + ") && " + ttjetsSingleLeptonCuts

def passTrigger(event, triggerNumList):
    """Checks if the event passed any trigger in the list"""
    passes = False
    for i in triggerNumList:
        if event.HLTDecision[i]:
            passes = True
            break
    return passes

def passSingleLeptonTrigger(event, isData=False, debug=False):
    if isData:
        trigNums = singleLeptonTriggerNumsData
    else:
        trigNums = singleLeptonTriggerNumsMC
    if debug: print("Event passes single lepton trigger")
    return passTrigger(event, trigNums)

def passDileptonTrigger(event, isData=False, debug=False):
    if debug: print("Note: dilepton trigger requirement is a pass-through right now")
    return True

def passHadronicTrigger(event, isData=False, debug=False):
    if debug: print("Note: hadronic trigger requirement is a pass-through right now")
    return True

def weight_mc(event, wHists, scale=1.0, debug=False):
    """Apply pileup weights and other known MC correction factors -- for razor control regions"""
    eventWeight = event.weight*scale
    #pileup reweighting
    if "pileup" in wHists:
        eventWeight *= wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(event.NPU_0))
    else:
        print("Error in standardRun2Weights: pileup reweighting histogram not found!")
        sys.exit()
    return eventWeight

def weight_data(event, wHists, scale=1.0, debug=False):
    eventWeight = scale
    if debug: print("Applying a weight of "+str(eventWeight))
    return eventWeight

def loadWeightHists(debug=False):
    """Returns a dict with necessary reweighting histograms"""
    puWeightFileName = "data/ScaleFactors/Placeholders/DummyRun2PileupWeights.root"
    muonWeightFileName = "data/ScaleFactors/Placeholders/DummyRun2MuonWeights.root"
    eleWeightFileName = "data/ScaleFactors/Placeholders/DummyRun2EleWeights.root"
    wHists = {}
    #pileup weight histogram
    if debug: print("Opening pileup weight file "+puWeightFileName)
    puFile = rt.TFile(puWeightFileName)
    wHists["pileup"] = puFile.Get("PUWeight_Run2")
    wHists["pileup"].SetDirectory(0)
    assert wHists["pileup"]
    #muon weight histogram
    if debug: print("Opening muon weight file "+muonWeightFileName)
    muonFile = rt.TFile(muonWeightFileName)
    wHists["muon"] = muonFile.Get("MuonWeight_Run2_Tight")
    wHists["muon"].SetDirectory(0)
    assert wHists["muon"]
    #ele weight histogram
    if debug: print("Opening ele weight file "+eleWeightFileName)
    eleFile = rt.TFile(eleWeightFileName)
    wHists["ele"] = eleFile.Get("EleWeight_Run2_Tight")
    wHists["ele"].SetDirectory(0)
    assert wHists["ele"]
    return wHists

