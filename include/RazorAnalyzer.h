//Class for analyzing ntuples produced by the RazorTuplizer framework
//
//Author: Caltech Razor team

#ifndef RazorAnalyzer_h
#define RazorAnalyzer_h

#include <RazorEvents.h> //This is a MakeClass of the RazorEvents tree in the ntuple to be analyzed

//ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TLorentzVector.h"

//C++ includes
#include <map>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class RazorAnalyzer: public RazorEvents {
    public :
        RazorAnalyzer(TTree *tree=0);
        virtual ~RazorAnalyzer();

        void EnableEventInfo();
        void EnableMuons();
        void EnableElectrons();
        void EnableTaus();
        void EnablePhotons();
        void EnableJets();
        void EnableFatJets();
        void EnableMet();
        void EnablePileup();
        void EnableMC();
        void EnableGenParticles();
        void EnableRazor();

        void EnableAll();

        //------ LIST OF ANALYSES ------//
        virtual void DummyAnalysis();
        virtual void RazorInclusive(string outFileName = "RazorInclusive.root", bool combineTrees = false);

        //functions in RazorAuxMuon.cc
	bool isVetoMuon(int i);
	bool isLooseMuon(int i);
        bool isTightMuon(int i);

        //functions in RazorAuxElectron.cc
        bool isVetoElectron(int i);
        bool isLooseElectron(int i);
        bool isTightElectron(int i);
        
        //functions in RazorAuxTau.cc
        bool isSelectedTau(int i);

        //functions in RazorAuxJet.cc
        bool isCSVL(int i);
        bool isCSVM(int i);
        bool isCSVT(int i);

        //functions in RazorAuxMisc.cc
	double deltaPhi(double phi1, double phi2);
	double deltaR(double eta1, double phi1, double eta2, double phi2);
        TLorentzVector makeTLorentzVector(double pt, double eta, double phi, double energy);
        vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
        double computeMR(TLorentzVector hem1, TLorentzVector hem2);
        double computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector met);
};

#endif
