//Class for analyzing ntuples produced by the RazorTuplizer framework
//
//Author: Caltech Razor team

#ifndef RazorAnalyzer_h
#define RazorAnalyzer_h

#include <RazorEvents.h> //This is a MakeClass of the RazorEvents tree in the ntuple to be analyzed
#include "FactorizedJetCorrector.h"

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
        void EnableIsoPFCandidates();
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
        virtual void HggRazor(string outFileName = "HggRazor.root", bool combineTrees = false);
        virtual void MatchedRazorInclusive(string outFileName = "MatchedRazorInclusive.root", bool combineTrees = false);
	virtual void RazorVetoLeptonStudy(string outputfilename = "RazorVetoLeptonStudy", bool combineTrees = false);
	virtual void ElectronNtupler(string outputfilename = "", int Option = -1);
	virtual void MuonNtupler(string outputfilename = "", int Option = -1);
	virtual void JetNtupler(string outputfilename = "", int Option = -1);
        virtual void RazorMetAna(string outFileName = "RazorMET.root");

        //functions in RazorAuxMuon.cc
	bool isVetoMuon(int i);
	bool isLooseMuon(int i);
        bool isTightMuon(int i);

        //functions in RazorAuxElectron.cc
        bool isVetoElectron(int i);
        bool isLooseElectron(int i);
        bool isTightElectron(int i);
	bool isMVANonTrigVetoElectron(int i);

        //functions in RazorAuxTau.cc
        bool isLooseTau(int i);
        bool isMediumTau(int i);
        bool isTightTau(int i);

        //functions in RazorAuxPhoton.cc
        bool photonPassesElectronVeto(int i);
        bool photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut);
        bool passesCutsBasedPhotonID(int i, double HoverECut, double SigmaIetaIetaCut, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut);
        bool passesCutsBasedPhotonIDNoIsoCuts(int i, double HoverECut, double SigmaIetaIetaCut);
        bool isLoosePhoton(int i);
        bool isMediumPhoton(int i);
        bool isMediumPhotonNoIsoCuts(int i);
        bool photonPassesMediumIsoCuts(int i);
        bool isTightPhoton(int i);

        //functions in RazorAuxJet.cc
	
        bool isOldCSVL(int i);
        bool isOldCSVM(int i);
        bool isOldCSVT(int i);
        bool isCSVL(int i);
        bool isCSVM(int i);
        bool isCSVT(int i);
	double JetEnergyCorrectionFactor( double jetRawPt, double jetEta, double jetPhi, double jetE,
					  double rho, double jetArea,
					  FactorizedJetCorrector *jetcorrector,  
					  bool printDebug = false);
	  
        //functions in RazorAuxMisc.cc
	double deltaPhi(double phi1, double phi2);
	double deltaR(double eta1, double phi1, double eta2, double phi2);
        TLorentzVector makeTLorentzVector(double pt, double eta, double phi, double energy);
	TLorentzVector makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass);
	vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
        double computeMR(TLorentzVector hem1, TLorentzVector hem2);
        double computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector met);

        bool passesHadronicRazorBaseline(double MR, double Rsq);
        bool passesLeptonicRazorBaseline(double MR, double Rsq);

        //enums
        enum RazorBox { //boxes for razor inclusive analysis
            MuEle, 
            MuMu,
            EleEle,
            MuMultiJet,
            MuJet,
            EleMultiJet,
            EleJet,
	    SoftLeptonMultiJet,
            MultiJet,
            TwoBJet,
            OneBJet,
            ZeroBJet,
            NONE
        };
};

#endif
