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
looseLeptonSignalRegionCuts = "MR > 400.000000 && Rsq > 0.250000 && mTLoose > 100 && nJets80 >= 2"
hadronicSignalRegionCuts = "MR > 400.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

razorCuts = {}
for box in razorBoxes:
    if box in hadronicRazorBoxes: 
        razorCuts[box] = appendBoxCuts(hadronicSignalRegionCuts, razorBoxes[box])
    elif box in looseLeptonRazorBoxes:
        razorCuts[box] = appendBoxCuts(looseLeptonSignalRegionCuts, razorBoxes[box])
    elif box in leptonicRazorBoxes:
        razorCuts[box] = appendBoxCuts(leptonicSignalRegionCuts, razorBoxes[box])
    elif box in dileptonRazorBoxes:
        razorCuts[box] = appendBoxCuts(dileptonSignalRegionCuts, razorBoxes[box])

leptonicSignalRegionBins = {
    "MR" : [300, 400, 500, 600, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }

hadronicSignalRegionBins = {
    "MR" : [400, 500, 600, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }

