import sys, os
import string
import numpy as np
import itertools
import ROOT as rt
import rootTools

##initialization

#check that the input file was specified
if len(sys.argv) < 2:
    print 'Usage: PhotonValidation.py filename.root'
    exit()

#switch ROOT to batch mode
rt.gROOT.SetBatch()
#turn on fit stats
rt.gStyle.SetOptFit(0111)

#get the path of this script and load the razor library
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname) #get the full path of this script
outpath = fullpath+"/output"
libpath = fullpath+"/lib"
rt.gSystem.Load(libpath+"/libRazorRun2.so")

print('Will compute photon selection efficiencies and validate the photon ID')

#load the TTree from the input file
inFile = rt.TFile(sys.argv[1])
events = inFile.Get("PhotonNtuple")

#plot photon efficiency as a function of pT, eta
c = rt.TCanvas("c", "c", 800, 600)
for level in ['Loose', 'Medium', 'Tight']:
    for ptCut in [0, 25, 40, 100]:
        #photon efficiency
        numerator = rt.TH1F("numerator", level+" Photon efficiency vs #eta, p_{T} > "+str(ptCut)+"; #eta; #epsilon", 50, -2.5, 2.5)
        denominator = rt.TH1F("denominator", level+" Photon efficiency vs #eta; #eta; #epsilon", 50, -2.5, 2.5)
        events.Draw("phoEta>>numerator", "phoMatchesGen && phoIs"+level+" && phoPt > "+str(ptCut))
        events.Draw("phoEta>>denominator", "phoMatchesGen && phoPt > "+str(ptCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"EfficiencyVsEtaForPt"+str(ptCut)+".pdf")
        #photon energy resolution
        c.SetLogy()
        EnergyResolution = rt.TH1F("EnergyResolution", level+" photon energy resolution in the barrel, p_{T} > "+str(ptCut)+"; #Delta E/E_{gen}", 50, 0., 1)
        events.Draw("deltaEOverEBest>>EnergyResolution", "phoMatchesGen && phoIs"+level+" && phoPt > "+str(ptCut))
        EnergyResolution.Draw()
        EnergyResolution.Fit("expo")
        c.Print(outpath+"/Photon"+level+"EnergyResolutionFoPt"+str(ptCut)+".pdf")
        c.SetLogy(rt.kFALSE)
        #photon isolation efficiency
        numerator = rt.TH1F("numerator", level+" Photon isolation efficiency vs #eta, p_{T} > "+str(ptCut)+"; #eta; #epsilon", 50, -2.5, 2.5)
        denominator = rt.TH1F("denominator", level+" Photon isolation efficiency vs #eta; #eta; #epsilon", 50, -2.5, 2.5)
        events.Draw("phoEta>>numerator", "phoMatchesGen && phoIsIsolated"+level+" && phoPt > "+str(ptCut))
        events.Draw("phoEta>>denominator", "phoMatchesGen && phoPt > "+str(ptCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"IsoEfficiencyVsEtaForPt"+str(ptCut)+".pdf")
    for etaCut in [2.5, 1.479]:
        #photon efficiency
        numerator = rt.TH1F("numerator", level+" Photon efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        denominator = rt.TH1F("denominator", level+" Photon efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        events.Draw("phoPt>>numerator", "phoMatchesGen && phoIs"+level+" && fabs(phoEta) < "+str(etaCut))
        events.Draw("phoPt>>denominator", "phoMatchesGen && fabs(phoEta) < "+str(etaCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"EfficiencyVsPtForEta"+str(etaCut)+".pdf")
        #photon isolation efficiency
        numerator = rt.TH1F("numerator", level+" Photon isolation efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        denominator = rt.TH1F("denominator", level+" Photon isolation efficiency vs p_{T}, |#eta| < "+str(etaCut)+"; p_{T}; #epsilon", 50, 0, 400)
        events.Draw("phoPt>>numerator", "phoMatchesGen && phoIsIsolated"+level+" && fabs(phoEta) < "+str(etaCut))
        events.Draw("phoPt>>denominator", "phoMatchesGen && fabs(phoEta) < "+str(etaCut))
        numerator.Divide(denominator)
        numerator.SetStats(0)
        numerator.SetMarkerStyle(20)
        numerator.Draw("pc")
        c.Print(outpath+"/Photon"+level+"IsoEfficiencyVsPtForEta"+str(etaCut)+".pdf")
