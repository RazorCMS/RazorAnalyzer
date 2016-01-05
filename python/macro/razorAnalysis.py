## Inclusive razor analysis definitions
import ROOT as rt
from array import array
import sys
import copy

#local imports
import macro

#####################################
### RAZOR CONTROL REGION DEFINITIONS
#####################################

#Hadronic trigger indices
hadronicTriggerNums = [134,135,136,137,138,139,140,141,142,143,144]

#Single lepton trigger indices
singleLeptonTriggerNumsData = [2,7,11,12,15,22,23,24,25,26,27,28,29]
singleLeptonTriggerNumsMC = [2,7,11,12,15,18,19,20,21,28,29,160]

#Dilepton trigger indices
dileptonTriggerNums = [41,43,30,31,47,48,49,50]

def appendTriggerCuts(cuts, trigNums):
    """Append a string of the form "(HLTDecision[t1] || HLTDecision[t2] || ... || HLTDecision[tN]) && " to the provided cut string, where t1...tN are the desired trigger numbers"""
    return '('+(' || '.join(['HLTDecision['+str(n)+']' for n in trigNums])) + ") && " + cuts

def appendBoxCuts(cuts, boxNums):
    """Append a string of the form "(box == b1 || box == b2 || ... || box == bN) && " to the provided cut string, where b1...bN are the desired box numbers"""
    return '('+(' || '.join(['box == '+str(n) for n in boxNums])) + ") && " + cuts

recommendedNoiseFilters = ["Flag_HBHENoiseFilter","Flag_HBHEIsoNoiseFilter","Flag_goodVertices","Flag_eeBadScFilter"]
def appendNoiseFilters(cuts, tree):
    ret = copy.copy(cuts)
    for bit in recommendedNoiseFilters:
        if hasattr(tree, bit):
            ret += " && " + bit
    return ret


### TTJets Single Lepton Control Region

#cuts
ttjetsSingleLeptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 25) || (abs(lep1Type) == 13 && lep1.Pt() > 20)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium > 0 && MR > 300 && Rsq > 0.15 && TMath::Finite(weight)"
ttjetsSingleLeptonCutsData = appendTriggerCuts(ttjetsSingleLeptonCuts, singleLeptonTriggerNumsData)
ttjetsSingleLeptonCutsMC = appendTriggerCuts(ttjetsSingleLeptonCuts, singleLeptonTriggerNumsMC)

#binning
ttjetsSingleLeptonBins = {
        "MR" : [300, 350, 400, 450, 500, 550, 700, 900, 4000],
        "Rsq": [0.15,0.175,0.20,0.225,0.25,0.30,0.41,1.5]
        }

### WJets Single Lepton Control Region

#cuts
wjetsSingleLeptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 25) || (abs(lep1Type) == 13 && lep1.Pt() > 20)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium == 0 && MR > 300 && Rsq > 0.15 && TMath::Finite(weight)"
wjetsSingleLeptonCutsData = appendTriggerCuts(wjetsSingleLeptonCuts, singleLeptonTriggerNumsData)
wjetsSingleLeptonCutsMC = appendTriggerCuts(wjetsSingleLeptonCuts, singleLeptonTriggerNumsMC)

#binning
wjetsSingleLeptonBins = {
        "MR" : [300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000],
        "Rsq": [0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5]
        }

### DYJets Dilepton Control Region

#cuts 
dyjetsDileptonCuts = "((abs(lep1Type) == 11 && abs(lep2Type) == 11) || (abs(lep1Type) == 13 && abs(lep2Type) == 13)) && lep1PassTight && lep2PassTight && lep1.Pt() > 30 && lep2.Pt() > 20 && mll > 80 && mll < 110 && NBJetsMedium == 0 && MR > 300 && Rsq > 0.15 && TMath::Finite(weight)"
dyjetsDileptonCutsData = appendTriggerCuts(dyjetsDileptonCuts, singleLeptonTriggerNumsData)
dyjetsDileptonCutsMC = appendTriggerCuts(dyjetsDileptonCuts, singleLeptonTriggerNumsMC)

#binning
dyjetsDileptonBins = {
        "MR" : [300, 350, 400, 450, 550],
        "Rsq": [0.15, 1.5]
        }

### TTJets Dilepton Control Region

#cuts 
ttjetsDileptonCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && (abs(lep2Type) == 11 || abs(lep2Type) == 13) && lep1PassTight && lep2PassTight && ((abs(lep1Type) == 11 && lep1.Pt() > 30) || (abs(lep1Type) == 13 && lep1.Pt() > 20)) && ((abs(lep2Type) == 11 && lep2.Pt() > 25) || (abs(lep2Type) == 13 && lep2.Pt() > 20)) && mll > 20 && ((abs(lep1Type) != abs(lep2Type)) || (mll < 76 || mll > 106)) && NBJetsMedium > 0 && MET > 40 && MR > 300 && Rsq > 0.15 && TMath::Finite(weight)"
ttjetsDileptonCutsData = appendTriggerCuts(ttjetsDileptonCuts, singleLeptonTriggerNumsData)
ttjetsDileptonCutsMC = appendTriggerCuts(ttjetsDileptonCuts, singleLeptonTriggerNumsMC)

#binning
ttjetsDileptonBins = {
    "MR" : [300, 350, 400, 450, 500, 550, 700],
    "Rsq": [0.15,0.175,0.20,0.225,0.25,0.30,1.5]
    }

### Veto Lepton Control Region
#vetoLeptonControlRegionCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 300 && Rsq > 0.15 && TMath::Finite(weight)"
vetoLeptonControlRegionCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && lep1PassVeto && lep1.Pt() > 5 && lep1MT > 30 && lep1MT < 100 && NJets80 >= 2 && MR > 300 && Rsq > 0.15 && TMath::Finite(weight)"
vetoLeptonControlRegionCutsData = appendTriggerCuts(vetoLeptonControlRegionCuts, hadronicTriggerNums)
vetoLeptonControlRegionCutsMC = appendTriggerCuts(vetoLeptonControlRegionCuts, hadronicTriggerNums)

vetoLeptonControlRegionBins = {
    "MR" : [400, 450, 500, 550, 700],
    "Rsq": [0.25,0.30,0.41,0.52, 0.64, 1.5]
    }


### WJets Single Lepton Invisible Control Region

#cuts
wjetsSingleLeptonInvCuts = "(abs(lep1Type) == 11 || abs(lep1Type) == 13) && ((abs(lep1Type) == 11 && lep1.Pt() > 25) || (abs(lep1Type) == 13 && lep1.Pt() > 20)) && MET > 30 && lep1MT > 30 && lep1MT < 100 && NBJetsMedium == 0 && NJets80_NoW >= 2 && MR_NoW > 300 && Rsq_NoW > 0.15 && TMath::Finite(weight)"
wjetsSingleLeptonInvCutsData = appendTriggerCuts(wjetsSingleLeptonInvCuts, singleLeptonTriggerNumsData)
wjetsSingleLeptonInvCutsMC = appendTriggerCuts(wjetsSingleLeptonInvCuts, singleLeptonTriggerNumsMC)

#binning
wjetsSingleLeptonInvBins = {
        "MR" : [300, 350, 400, 450, 500, 550, 700, 900, 1200, 4000],
        "Rsq": [0.15,0.175,0.20,0.225, 0.25,0.30,0.41,0.52,1.5]
        }

### DYJets Dilepton Invisible Control Region

#cuts 
dyjetsDileptonInvCuts = "((abs(lep1Type) == 11 && abs(lep2Type) == 11) || (abs(lep1Type) == 13 && abs(lep2Type) == 13)) && lep1.Pt() > 30 && lep2.Pt() > 20 && mll > 80 && mll < 110 && NBJetsMedium == 0 && NJets80_NoZ >= 2 && MR_NoZ > 300 && Rsq_NoZ > 0.15 && TMath::Finite(weight)"
dyjetsDileptonInvCutsData = appendTriggerCuts(dyjetsDileptonInvCuts, singleLeptonTriggerNumsData)
dyjetsDileptonInvCutsMC = appendTriggerCuts(dyjetsDileptonInvCuts, singleLeptonTriggerNumsMC)

#binning
dyjetsDileptonInvBins = {
        "MR" : [300, 350, 400, 450, 550],
        "Rsq": [0.15, 1.5]
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
dileptonSignalRegionCuts = "MR > 400.000000 && Rsq > 0.150000 && abs(dPhiRazor) < 2.8"
leptonicSignalRegionCuts = "MR > 400.000000 && Rsq > 0.150000 && mT > 120"
looseLeptonSignalRegionCuts = "MR > 500.000000 && Rsq > 0.250000 && mTLoose > 100 && nJets80 >= 2"
hadronicSignalRegionCuts = "MR > 500.000000 && Rsq > 0.250000 && abs(dPhiRazor) < 2.8 && nJets80 >= 2"

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

#NOTE: razor binning should generally be taken directly from a config, to avoid conflicts with fit code
leptonicSignalRegionBins = {
    "MR" : [300, 400, 500, 600, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }
leptonicSidebandBins = {
    "MR" : [300, 350, 450, 550, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.15,0.20,0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }
leptonicBlindBins = [(x,y) for x in range(3,len(leptonicSidebandBins["MR"])+1) for y in range(2,len(leptonicSidebandBins["Rsq"])+1)]
hadronicSignalRegionBins = {
    "MR" : [400, 500, 600, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }
hadronicSidebandBins = {
    "MR" : [400, 450, 550, 700, 900, 1200, 1600, 2500, 4000],
    "Rsq" : [0.25,0.30,0.41,0.52,0.64,0.8,1.5]
    }
hadronicBlindBins = [(x,y) for x in range(3,len(hadronicSidebandBins["MR"])+1) for y in range(2,len(hadronicSidebandBins["Rsq"])+1)]
