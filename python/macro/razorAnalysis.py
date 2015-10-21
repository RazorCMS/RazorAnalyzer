## Inclusive razor analysis definitions
import ROOT as rt
from array import array
import sys

#local imports
import macro

#####################################
### RAZOR CONTROL REGION DEFINITIONS
#####################################

#Hadronic trigger indices
hadronicTriggerNums = [134,135,136,137,138,139,140,141,142,143,144]

#Single lepton trigger indices
singleLeptonTriggerNumsData = [2,7,11,12,15,22,23,24,25,26,27,28,29]
singleLeptonTriggerNumsMC = [2,7,11,12,15,18,19,20,21,28,29]

#Dilepton trigger indices
dileptonTriggerNums = [41,43,30,31,47,48,49,50]

def appendTriggerCuts(cuts, trigNums):
    """Append a string of the form "(HLTDecision[t1] || HLTDecision[t2] || ... || HLTDecision[tN]) && " to the provided cut string, where t1...tN are the desired trigger numbers"""
    return '('+(' || '.join(['HLTDecision['+str(n)+']' for n in trigNums])) + ") && " + cuts

def appendBoxCuts(cuts, boxNums):
    """Append a string of the form "(box == b1 || box == b2 || ... || box == bN) && " to the provided cut string, where b1...bN are the desired box numbers"""
    return '('+(' || '.join(['box == '+str(n) for n in boxNums])) + ") && " + cuts

### TTJets Single Lepton Control Region

#cuts
ttjetsSingleLeptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1.Pt() > 30 && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium > 0 && MR > 300 && Rsq > 0.15"
ttjetsSingleLeptonCutsData = appendTriggerCuts(ttjetsSingleLeptonCuts, singleLeptonTriggerNumsData)
ttjetsSingleLeptonCutsMC = appendTriggerCuts(ttjetsSingleLeptonCuts, singleLeptonTriggerNumsMC)

#binning
ttjetsSingleLeptonBins = {
        "MR" : [300, 350, 400, 450, 500, 550, 700, 900, 4000],
        "Rsq": [0.15,0.175,0.20,0.225,0.25,0.30,0.41, 1.5]
        }

### WJets Single Lepton Control Region

#cuts
wjetsSingleLeptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && lep1.Pt() > 30 && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium == 0 && MR > 300 && Rsq > 0.15"
wjetsSingleLeptonCutsData = appendTriggerCuts(wjetsSingleLeptonCuts, singleLeptonTriggerNumsData)
wjetsSingleLeptonCutsMC = appendTriggerCuts(wjetsSingleLeptonCuts, singleLeptonTriggerNumsMC)

#binning
wjetsSingleLeptonBins = {
        "MR" : [300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000],
        "Rsq": [0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5]
        }

### DYJets Dilepton Control Region

#cuts (NOTE: requiring tight leptons at the moment)
dyjetsDileptonCuts = "((abs(lep1Type) == 11 && abs(lep2Type) == 11) || (abs(lep1Type) == 13 && abs(lep2Type) == 13)) && lep1PassLoose && lep2PassLoose && lep1.Pt() > 25 && lep2.Pt() > 25 && mll > 80 && mll < 110 && NBJetsMedium == 0 && MR > 300 && Rsq > 0.15"
dyjetsDileptonCutsData = appendTriggerCuts(dyjetsDileptonCuts, singleLeptonTriggerNumsData)
dyjetsDileptonCutsMC = appendTriggerCuts(dyjetsDileptonCuts, singleLeptonTriggerNumsMC)

#binning
dyjetsDileptonBins = {
        "MR" : [300, 350, 400, 450, 550],
        "Rsq": [0.15, 1.5]
        }

### TTJets Dilepton Control Region

#cuts (NOTE: requiring tight leptons at the moment)
ttjetsDileptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 30 && mll > 20 && ((abs(lep1Type) != abs(lep2Type)) || (mll > 76 && mll < 106)) && NBJetsMedium > 0 && MET > 40 && MR > 300 && Rsq > 0.15"
ttjetsDileptonCutsData = appendTriggerCuts(ttjetsDileptonCuts, singleLeptonTriggerNumsData)
ttjetsDileptonCutsMC = appendTriggerCuts(ttjetsDileptonCuts, singleLeptonTriggerNumsMC)

#binning
ttjetsDileptonBins = {
    "MR" : [300, 350, 400, 450, 500, 550, 700],
    "Rsq": [0.15,0.175,0.20,0.225,0.25,0.30,1.5]
    }

### Signal region

#boxes
razorBoxes = {
        "MuEle" : [0], 
        "MuMu" : [1],
        "EleEle" : [2],
        "MuSixJet" : [3],
        "MuFourJet" : [4],
        "MuJet" : [5],
        "EleSixJet" : [6],
        "EleFourJet" : [7],
        "EleJet" : [8],
        "LooseLeptonSixJet" : [9],
        "LooseLeptonFourJet" : [10],
        "LooseLeptonDiJet" : [13],
        "SixJet" : [11],
        "FourJet" : [12],
        "DiJet" : [14],	  
        "MuMultiJet" : [3,4,18],
        "EleMultiJet" : [6,7,19],
        "LooseLeptonMultiJet" : [9,10,20],
        "MultiJet" : [11,12,21],
        }
hadronicRazorBoxes = ["DiJet", "FourJet", "SixJet", "MultiJet"]
looseLeptonRazorBoxes = ["LooseLeptonDiJet", "LooseLeptonFourJet", "LooseLeptonSixJet", "LooseLeptonMultiJet"]
leptonicRazorBoxes = ["MuJet", "MuFourJet", "MuSixJet", "MuMultiJet","EleJet", "EleFourJet", "EleSixJet", "EleMultiJet"]
dileptonRazorBoxes = ["MuEle", "MuMu", "EleEle"]

#cuts 
dileptonSignalRegionCuts = "MR > 300.000000 && Rsq > 0.150000 && abs(dPhiRazor) < 2.8"
leptonicSignalRegionCuts = "MR > 300.000000 && Rsq > 0.150000 && mT > 100"
looseLeptonSignalRegionCuts = "MR > 400.000000 && Rsq > 0.250000 && mT > 100"
#looseLeptonSignalRegionCuts = "MR > 400.000000 && Rsq > 0.250000 && mT > 100 && nJets80 >= 2"
hadronicSignalRegionCuts = "MR > 400.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8"
#hadronicSignalRegionCuts = "MR > 400.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

razorCutsData = {}
razorCutsMC = {}
for box in razorBoxes:
    if box in hadronicRazorBoxes: 
        tmp = appendBoxCuts(hadronicSignalRegionCuts, razorBoxes[box])
        razorCutsData[box] = appendTriggerCuts(tmp, hadronicTriggerNums)
        razorCutsMC[box] = appendTriggerCuts(tmp, hadronicTriggerNums)
    elif box in looseLeptonRazorBoxes:
        tmp = appendBoxCuts(looseLeptonSignalRegionCuts, razorBoxes[box])
        razorCutsData[box] = appendTriggerCuts(tmp, hadronicTriggerNums)
        razorCutsMC[box] = appendTriggerCuts(tmp, hadronicTriggerNums)
    elif box in leptonicRazorBoxes:
        tmp = appendBoxCuts(leptonicSignalRegionCuts, razorBoxes[box])
        razorCutsData[box] = appendTriggerCuts(tmp, singleLeptonTriggerNumsData)
        razorCutsMC[box] = appendTriggerCuts(tmp, singleLeptonTriggerNumsMC)
    elif box in dileptonRazorBoxes:
        tmp = appendBoxCuts(dileptonSignalRegionCuts, razorBoxes[box])
        razorCutsData[box] = appendTriggerCuts(tmp, dileptonTriggerNums)
        razorCutsMC[box] = appendTriggerCuts(tmp, dileptonTriggerNums)

leptonicSignalRegionBins = {
    "MR" : [300, 400, 500, 600, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }

hadronicSignalRegionBins = {
    "MR" : [400, 500, 600, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }

#####################################
### WEIGHT AND TRIGGER INFO
#####################################

def passTrigger(event, triggerNumList):
    """Checks if the event passed any trigger in the list"""
    passes = False
    for i in triggerNumList:
        if event.HLTDecision[i]:
            passes = True
            break
    return passes

def passSingleLeptonTrigger(event, isData=False, debugLevel=0):
    if isData:
        trigNums = singleLeptonTriggerNumsData
    else:
        trigNums = singleLeptonTriggerNumsMC
    if debugLevel > 1: print("Event passes single lepton trigger")
    return passTrigger(event, trigNums)

def passDileptonTrigger(event, isData=False, debugLevel=0):
    if debugLevel > 1: print("Note: dilepton trigger requirement is a pass-through right now")
    return True

def passHadronicTrigger(event, isData=False, debugLevel=0):
    if debugLevel > 1: print("Note: hadronic trigger requirement is a pass-through right now")
    return True

def weight_mc_signalregion(event, wHists, scale=1.0, opts=[], debugLevel=0):
    """Assume pileup weight, lepton efficiency weight, and trigger weight, have been applied"""
    eventWeight = event.weight*scale
    if debugLevel > 1: 
        print "Weight from ntuple:",event.weight
        print "Scale by:",scale

    if debugLevel > 1: 
        print "event weight:",eventWeight
    return eventWeight

def weight_mc(event, wHists, scale=1.0, opts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], debugLevel=0):
    """Apply pileup weights and other known MC correction factors -- for razor control regions"""
    eventWeight = event.weight*scale
    if debugLevel > 1: 
        print "Weight from ntuple:",event.weight
        print "Scale by:",scale
    #pileup reweighting
    if "doPileupWeights" in opts:
        if "pileup" in wHists:
            pileupWeight = wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(event.NPV))
            #pileupWeight = wHists["pileup"].GetBinContent(wHists["pileup"].GetXaxis().FindFixBin(event.NPU_0))
            if debugLevel > 1: print "NPV weight:",pileupWeight,"( NPV=",event.NPV,")"
            eventWeight *= pileupWeight
        else:
            print "Error in weight_mc: pileup reweighting histogram not found!"
            sys.exit()

    #lepton scale factors
    if "doLep1Weights" in opts:
        if not ("ele" in wHists):
            print "Error in weight_mc: electron scale factor histogram not found!"
            sys.exit()
        if not ("muon" in wHists):
            print "Error in weight_mc: muon scale factor histogram not found!"
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
        if "doLep2Weights" in opts:
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
        eventWeight *= leptonWeight

    #trigger scale factors
    if "do1LepTrigWeights" in opts:
        if not ("eletrig" in wHists):
            print "Error in weight_mc: electron trigger scale factor histogram not found!"
            sys.exit()
        if not ("muontrig" in wHists):
            print "Error in weight_mc: muon trigger scale factor histogram not found!"
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
        if "doLep2TrigWeights" in opts:
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
            print "1-lepton trigger weight:",trigWeight
        eventWeight *= trigWeight

    if debugLevel > 1: 
        print "event weight:",eventWeight
    return eventWeight

def weight_data(event, wHists, scale=1.0, opts=[], debugLevel=0):
    eventWeight = scale
    if debugLevel > 1: print("Applying a weight of "+str(eventWeight))
    return eventWeight

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

def loadScaleFactorHists(sfFilename="RazorScaleFactors.root", processNames=[], debugLevel=0):
    """Returns a dict with available scale factor histograms"""
    sfFile = rt.TFile.Open(sfFilename)
    sfHists = {}
    for name in processNames:
        if debugLevel > 0: print "Looking for scale factor histogram for",name,"...",
        tmp = sfFile.Get(name+"ScaleFactors")
        if tmp: 
            sfHists[name] = tmp
            sfHists[name].SetDirectory(0)
            if debugLevel > 0: print "Found!"
        else:
            if debugLevel > 0: print ""
    sfFile.Close()
    return sfHists

#####################################
### PLOTTING/SCALE FACTOR MACROS
#####################################

titles = {
    "MR": "M_{R} (GeV)", 
    "Rsq": "R^{2}",
    "mll": "m_{ll} (GeV)",
    }
def makeControlSampleHists(regionName="TTJetsSingleLepton", filenames={}, samples=[], cutsMC="", cutsData="", bins={}, logX=True, lumiMC=1, lumiData=3000, weightHists={}, sfHists={}, treeName="ControlSampleEvent",dataName="Data", opts=["doPileupWeights", "doLep1Weights", "do1LepTrigWeights"], debugLevel=0):
    #setup files and trees
    inputs = filenames
    files = {name:rt.TFile.Open(inputs[name]) for name in inputs} #get input files
    for name in inputs: 
        assert files[name] #check input files
        if debugLevel > 0: print "Opened file",inputs[name]
    trees = macro.makeTreeDict(files, treeName, debugLevel)

    #define histograms to fill
    hists = {name:{} for name in inputs}
    for name in inputs:
        #1D histograms
        for var in bins:
            if var in titles: title=titles[var]
            else: title = var
            hists[name][var] = rt.TH1F(regionName+var+name, title, len(bins[var])-1, array('d',bins[var]))
        #2D MR-Rsq histogram
        if "MR" in bins and "Rsq" in bins:
            hists[name][("MR","Rsq")] = rt.TH2F(regionName+"MRRsq"+name, "R^{2} vs M_{R}", len(bins["MR"])-1, array('d',bins["MR"]), len(bins["Rsq"])-1, array('d',bins["Rsq"]))
        for var in hists[name]: 
            hists[name][var].Sumw2()
            hists[name][var].SetDirectory(0)
    listOfVars = hists.itervalues().next().keys() #list of the variable names
    
    #fill histograms
    print("Data:")
    macro.loopTree(trees[dataName], weightF=weight_data, cuts=cutsData, hists=hists[dataName], weightHists=weightHists, opts=[], debugLevel=debugLevel) 

    print("MC:")
    macro.loopTrees(trees, weightF=weight_mc, cuts=cutsMC, hists={name:hists[name] for name in samples}, weightHists=weightHists, sfHists=sfHists, scale=lumiData*1.0/lumiMC, opts=opts, debugLevel=debugLevel) 

    #print histograms
    c = rt.TCanvas(regionName+"c", regionName+"c", 800, 600)
    rt.SetOwnership(c, False)
    macro.basicPrint(hists, mcNames=samples, varList=listOfVars, c=c, printName=regionName, logx=logX, dataName=dataName, ymin=0.1, lumistr=str(lumiData)+" pb^{-1}")

    #close files and return
    for f in files: files[f].Close()
    return hists

def appendScaleFactors(process="TTJets", hists={}, sfHists={}, var=("MR","Rsq"), dataName="Data", normErrFraction=0.2, printTable=True, debugLevel=0):
    """Subtract backgrounds and make the data/MC histogram for the given process.
    Also makes up/down histograms corresponding to 20% shifts in the background normalization"""
    if debugLevel > 0: 
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
    bgProcesses = [mcProcess for mcProcess in hists if mcProcess != process and mcProcess != dataName]
    for mcProcess in bgProcesses:
        if var not in hists[mcProcess]:
            print "Error in appendScaleFactors: could not find",var," in hists[",mcProcess,"]!"
            return
        if debugLevel > 0: print "Subtracting",mcProcess,"from",dataName,"distribution"
        sfHists[process].Add(hists[mcProcess][var], -1) 
        if mcProcess not in sfHists: #background normalization uncertainty affects only processes with no scale factors
            if debugLevel > 0: print "Process",mcProcess,"has no associated scale factors.  Its normalization will be given a 20% uncertainty"
            sfHists[process+"NormUp"].Add(hists[mcProcess][var], -(1+normErrFraction)) 
            sfHists[process+"NormDown"].Add(hists[mcProcess][var], -(1/(1+normErrFraction))) 
        else: 
            if debugLevel > 0: print "Process",mcProcess,"has associated scale factors.  No further uncertainty will be applied to its normalization."
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

    if debugLevel > 0:
        print "Scale factor histograms after adding",process,":"
        print sfHists

    if printTable:
        xbinLowEdges = []
        xbinUpEdges = []
        ybinLowEdges = []
        ybinUpEdges = []
        sysUncerts = {mcProcess:[] for mcProcess in bgProcesses}
        statUncerts = []
        sfs = []
        #for each bin, get values for all table columns
        for bx in range(1, sfHists[process].GetNbinsX()+1):
            for by in range(1, sfHists[process].GetNbinsY()+1):
                xbinLowEdges.append('%.0f' % (sfHists[process].GetXaxis().GetBinLowEdge(bx)))
                xbinUpEdges.append('%.0f' % (sfHists[process].GetXaxis().GetBinUpEdge(bx)))
                ybinLowEdges.append(str(sfHists[process].GetYaxis().GetBinLowEdge(by)))
                ybinUpEdges.append(str(sfHists[process].GetYaxis().GetBinUpEdge(by)))
                scaleFactor = sfHists[process].GetBinContent(bx,by)
                sfs.append('%.3f' % (scaleFactor))
                if scaleFactor > 0:
                    statUncerts.append('%.1f\\%%' % (100*(sfHists[process].GetBinError(bx,by)/scaleFactor)))
                else: 
                    statUncerts.append('--')
                for mcProcess in bgProcesses: 
                    dataYield = hists[dataName][var].GetBinContent(bx,by)
                    if dataYield > 0 and scaleFactor > 0:
                        sysUncerts[mcProcess].append('%.1f\\%%' % (100*abs(hists[mcProcess][var].GetBinContent(bx,by)*normErrFraction*1.0/dataYield/scaleFactor)))
                    else: 
                        sysUncerts[mcProcess].append('--')
        xRanges = [low+'-'+high for (low, high) in zip(xbinLowEdges, xbinUpEdges)]
        yRanges = [low+'-'+high for (low, high) in zip(ybinLowEdges, ybinUpEdges)]
        headers=[var[0], var[1], "Scale factor", "Stat.\\ unc."]
        cols = [xRanges, yRanges, sfs, statUncerts]
        for mcProcess in bgProcesses: 
            headers.extend(["Unc.\\ from "+mcProcess])
            cols.extend([sysUncerts[mcProcess]])
        macro.table_basic(headers, cols, caption="Scale factors for "+process+" background", printstr="scaleFactorTable"+process)
