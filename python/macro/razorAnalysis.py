## Inclusive razor analysis definitions
import ROOT as rt
from array import array

#local imports
import macro

#Single lepton trigger indices
#singleLeptonTriggerNumsData = [1,2,8,20,22,24,25] #used in razor inclusive
#singleLeptonTriggerNumsMC = [1,2,8,17,18,19,24,25] #used in razor inclusive
singleLeptonTriggerNumsData = [3,8,11,12,14,21,22,23,24,25,26,27]
singleLeptonTriggerNumsMC = [3,8,11,12,14,17,18,19,20,26,27]

#Dilepton trigger indices (NOTE: using single lepton paths for now)
dileptonTriggerNumsData = singleLeptonTriggerNumsData
dileptonTriggerNumsMC = singleLeptonTriggerNumsMC

### TTJets Single Lepton Control Region

#cuts
ttjetsSingleLeptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1Pt > 30 && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium > 0 && MR > 300 && Rsq > 0.15"
ttjetsSingleLeptonCutsData = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in singleLeptonTriggerNumsData])) + ") && " + ttjetsSingleLeptonCuts
ttjetsSingleLeptonCutsMC = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in singleLeptonTriggerNumsMC])) + ") && " + ttjetsSingleLeptonCuts

#binning
ttjetsSingleLeptonMRBins = [300, 350, 400, 450, 500, 550, 700, 900, 4000]
ttjetsSingleLeptonRsqBins = [0.15,0.175,0.20,0.225,0.25,0.30,0.41, 1.5]

### WJets Single Lepton Control Region

#cuts
wjetsSingleLeptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1Pt > 30 && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium == 0 && MR > 300 && Rsq > 0.15"
wjetsSingleLeptonCutsData = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in singleLeptonTriggerNumsData])) + ") && " + wjetsSingleLeptonCuts
wjetsSingleLeptonCutsMC = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in singleLeptonTriggerNumsMC])) + ") && " + wjetsSingleLeptonCuts

#binning
wjetsSingleLeptonMRBins = [300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000]
wjetsSingleLeptonRsqBins = [0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5]

### DYJets Dilepton Control Region

#cuts (NOTE: requiring tight leptons at the moment)
dyjetsDileptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1PassLoose && lep2PassLoose && lep1.Pt() > 25 && lep2.Pt() > 25 && abs(lep1Type) == abs(lep2Type) && mll > 80 && mll < 110 && NBJetsMedium == 0 && MR > 300 && Rsq > 0.15"
dyjetsDileptonCutsData = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in dileptonTriggerNumsData])) + ") && " + dyjetsDileptonCuts
dyjetsDileptonCutsMC = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in dileptonTriggerNumsMC])) + ") && " + dyjetsDileptonCuts

#binning
dyjetsDileptonMRBins = [300, 350, 400, 450, 550]
dyjetsDileptonRsqBins = [0.15, 1.5]

### TTJets Dilepton Control Region

#cuts (NOTE: requiring tight leptons at the moment)
ttjetsDileptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 30 && mll > 20 && ((abs(lep1Type) != abs(lep2Type)) || (mll > 76 && mll < 106)) && NBJetsMedium > 0 && MET > 40 && MR > 300 && Rsq > 0.15"
ttjetsDileptonCutsData = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in dileptonTriggerNumsData])) + ") && " + ttjetsDileptonCuts
ttjetsDileptonCutsMC = '('+(' || '.join(['HLTDecision['+str(n)+']' for n in dileptonTriggerNumsMC])) + ") && " + ttjetsDileptonCuts

#binning
ttjetsDileptonMRBins = [300, 350, 400, 450, 500, 550, 700]
ttjetsDileptonRsqBins = [0.15,0.175,0.20,0.225,0.25,0.30,1.5]

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
    #if debug: print("Applying a weight of "+str(eventWeight))
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

    puFile.Close()
    muonFile.Close()
    eleFile.Close()
    return wHists

def makeControlSampleHists(regionName="TTJetsSingleLepton", filenames={}, samples=[], cutsMC="", cutsData="", binsMR=[400,4000], binsRsq=[0.15,1.5], lumiMC=1, lumiData=3000, weightHists={}, sfHists={}, treeName="ControlSampleEvent", debug=False):
    #setup files and trees
    inputs = filenames
    files = {name:rt.TFile.Open(inputs[name]) for name in inputs} #get input files
    for name in inputs: assert files[name] #check input files
    trees = macro.makeTreeDict(files, treeName, debug)

    #define histograms to fill
    hists = {name:{} for name in inputs}
    for name in inputs:
        hists[name]["MR"] = rt.TH1F(regionName+"MR"+name, "M_{R} (GeV)", len(binsMR)-1, array('d',binsMR))
        hists[name]["Rsq"] = rt.TH1F(regionName+"Rsq"+name, "R^{2}", len(binsRsq)-1, array('d',binsRsq))
        hists[name][("MR","Rsq")] = rt.TH2F(regionName+"MRRsq"+name, "R^{2} vs M_{R}", len(binsMR)-1, array('d',binsMR), len(binsRsq)-1, array('d',binsRsq))
        for var in hists[name]: 
            hists[name][var].Sumw2()
            hists[name][var].SetDirectory(0)
    listOfVars = hists.itervalues().next().keys() #list of the variable names
    
    #fill histograms
    print("MC:")
    macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:hists[name] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, debug=debug) 

    print("Data:")
    macro.loopTree(trees["Data"], weightF=weight_data, cuts=cutsData, hists=hists["Data"], weightHists=weightHists, debug=debug) 

    #print histograms
    c = rt.TCanvas(regionName+"c", regionName+"c", 800, 600)
    rt.SetOwnership(c, False)
    macro.basicPrint(hists, mcNames=samples, varList=listOfVars, c=c, printName=regionName, logx=True)

    #close files and return
    for f in files: files[f].Close()
    return hists

def appendScaleFactors(process="TTJets", hists={}, sfHists={}, var=("MR","Rsq"), dataName="Data", normErrFraction=0.2, debug=False):
    """Subtract backgrounds and make the data/MC histogram for the given process.
    Also makes up/down histograms corresponding to 20% shifts in the background normalization"""
    if debug: 
        print "Scale factor histograms so far:"
        print sfHists
    #warn if this scale factor histogram already exists
    if process in sfHists:
        print "Warning in appendScaleFactors: ",process," scale factor histogram already exists!  Will overwrite..."
    #get the target MC histogram
    if process not in hists:
        print "Error in appendScaleFactors: target MC histogram (",process,") was not found!"
        return
    if var not in hists[process]:
        print "Error in appendScaleFactors: could not find ",var," in hists[",process,"]!"
        return
    if dataName not in hists:
        print "Error in appendScaleFactors: target data histogram (",dataName,") was not found!"
        return
    if var not in hists[dataName]:
        print "Error in appendScaleFactors: could not find ",var," in hists[",dataName,"]!"
    print "Making scale factor histogram for ",process
    sfHists[process] = hists[dataName][var].Clone(process+"ScaleFactors")
    sfHists[process].SetDirectory(0)
    #make systematic error histograms
    sfHists[process+"NormUp"] = hists[dataName][var].Clone(process+"ScaleFactorsUp")
    sfHists[process+"NormUp"].SetDirectory(0)
    sfHists[process+"NormDown"] = hists[dataName][var].Clone(process+"ScaleFactorsDown")
    sfHists[process+"NormDown"].SetDirectory(0)

    #subtract backgrounds in data
    for mcProcess in hists:
        if mcProcess == process: continue
        if mcProcess == dataName: continue
        if var not in hists[mcProcess]:
            print "Error in appendScaleFactors: could not find",var," in hists[",mcProcess,"]!"
            return
        if debug: print "Subtracting",mcProcess,"from",dataName,"distribution"
        sfHists[process].Add(hists[mcProcess][var], -1) 
        if mcProcess not in sfHists: #background normalization uncertainty affects only processes with no scale factors
            if debug: print "Process",mcProcess,"has no associated scale factors.  Its normalization will be given a 20% uncertainty"
            sfHists[process+"NormUp"].Add(hists[mcProcess][var], -(1+normErrFraction)) 
            sfHists[process+"NormDown"].Add(hists[mcProcess][var], -(1/(1+normErrFraction))) 
        else: 
            if debug: print "Process",mcProcess,"has associated scale factors.  No further uncertainty will be applied to its normalization."
            sfHists[process+"NormUp"].Add(hists[mcProcess][var], -1) 
            sfHists[process+"NormDown"].Add(hists[mcProcess][var], -1) 

    #divide data/MC
    sfHists[process].Divide(hists[process][var])
    sfHists[process+"NormUp"].Divide(hists[process][var])
    sfHists[process+"NormDown"].Divide(hists[process][var])

    #zero any negative scale factors
    for bx in range(1, sfHists[process].GetNbinsX()+1):
        for by in range(1, sfHists[process].GetNbinsY()+1):
            sfHists[process].SetBinContent(bx,by,max(0., sfHists[process].GetBinContent(bx,by)))
            sfHists[process+"NormUp"].SetBinContent(bx,by,max(0., sfHists[process+"NormUp"].GetBinContent(bx,by)))
            sfHists[process+"NormDown"].SetBinContent(bx,by,max(0., sfHists[process+"NormDown"].GetBinContent(bx,by)))

    if debug:
        print "Scale factor histograms after adding",process,":"
        print sfHists
