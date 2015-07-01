//Class for analyzing ntuples produced by the RazorTuplizer framework
//
//Author: Caltech Razor team

#ifndef RazorAnalyzer_h
#define RazorAnalyzer_h

#include "RazorEvents.h" //This is a MakeClass of the RazorEvents tree in the ntuple to be analyzed
#include "FactorizedJetCorrector.h"
#include "SimpleJetResolution.h"
#include "HggRazorAuxPhoton.hh"

//ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include "TLorentzVector.h"
#include "TRandom3.h"

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
        virtual void RazorInclusive(string outFileName = "RazorInclusive.root", bool combineTrees = false, bool isData = false, bool isRunOne = false);
        virtual void HggRazor(string outFileName = "HggRazor.root", bool combineTrees = false);//Run1
	virtual void HggRazorRun2(string outFileName = "HggRazorRun2.root", bool combineTrees = false);//Run2
        virtual void MatchedRazorInclusive(string outFileName = "MatchedRazorInclusive.root", bool combineTrees = false);
	virtual void RazorVetoLeptonStudy(string outputfilename = "RazorVetoLeptonStudy", bool combineTrees = false);
	virtual void ElectronNtupler(string outputfilename = "", int Option = -1);
	virtual void MuonNtupler(string outputfilename = "", int Option = -1);
	virtual void TauNtupler(string outputfilename = "", int Option = -1);
	virtual void JetNtupler(string outputfilename = "", int Option = -1);
        virtual void PhotonNtupler(string outputfilename = "PhotonNtuple.root", int Option = -1); 
        virtual void RazorMetAna(string outFileName = "RazorMET.root");
	virtual void RazorDM(string outFileName = "RazorInclusive.root", bool combineTrees = false);
	virtual void RazorControlRegions(string outFileName = "RazorControlRegions.root", int option = -1, bool isData = false, bool isRunOne = false);
	virtual void VetoLeptonEfficiencyControlRegion(string outFileName = "TTBarTagAndProbeRegion.root", int option = 0);
        virtual void RazorPhotonStudy(string outputfilename = "RazorPhotonStudy.root", bool isData = false, bool filterEvents = true, bool isRunOne = true);
        virtual void MakeMCPileupDistribution(string outputfilename = "MCPileupDistribution.root", string label = "defaultSample");
	virtual void RazorZAnalysis(string outFileName = "RazorZAnalysis.root", bool combineTrees = false);
        virtual void HbbRazor(string outFileName = "HbbRazor.root", bool combineTrees = false, bool isData = false, bool isRunOne = false);
        virtual void HZZRazor(string outFileName = "HZZRazor.root", bool isData = false, bool isRunOne = false);

        //functions in RazorAuxMuon.cc
	bool isVetoMuon(int i);
	bool isLooseMuon(int i);
        bool isTightMuon(int i);
	bool passVetoMuonID(int i);
	bool passLooseMuonID(int i);
        bool passTightMuonID(int i);
	bool passVetoMuonIso(int i);
	bool passLooseMuonIso(int i);
        bool passTightMuonIso(int i);
        bool passRunOneHZZMuonPreselection(int i);
        bool isRunOneHZZMuon(int i);

        //functions in RazorAuxElectron.cc
        bool isVetoElectron(int i);
        bool isLooseElectron(int i);
        bool isTightElectron(int i);
	bool isMVANonTrigVetoElectron(int i);
        bool passTightElectronID(int i);
        bool passLooseElectronID(int i);
	bool passMVANonTrigVetoElectronID(int i);
        bool passTightElectronIso(int i);
        bool passLooseElectronIso(int i);
	bool passMVANonTrigVetoElectronIso(int i);
	bool passRunOneHZZElectronIso(int i);
	bool passRunOneHZZElectronPreselection(int i);
	bool isRunOneHZZElectron(int i);
        bool isRunOneLooseElectron(int i);
        bool isRunOneTightElectron(int i);
	double getElectronEfficiencyRunOne( string selectionType, double pt, double eta);

        //functions in RazorAuxTau.cc
        bool isLooseTau(int i);
        bool isMediumTau(int i);
        bool isTightTau(int i);

        //functions in RazorAuxPhoton.cc
        bool photonPassesElectronVeto(int i);
        bool photonPassesIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut);
        bool photonPassesRunOneIsolation(int i, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut);
        bool passesCutsBasedPhotonID(int i, double HoverECut, double SigmaIetaIetaCut, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut);
        bool passesRunOneCutsBasedPhotonID(int i, double HoverECut, double SigmaIetaIetaCut, double PFChHadIsoCut, double PFNeuHadIsoCut, double PFPhotIsoCut);
        bool passesCutsBasedPhotonIDNoIsoCuts(int i, double HoverECut, double SigmaIetaIetaCut);
        bool isLoosePhoton(int i);
        bool isMediumPhoton(int i);
        bool isMediumRunOnePhoton(int i);
        bool isTightRunOnePhoton(int i);
        bool isMediumPhotonNoIsoCuts(int i);
        bool photonPassesLooseIsoCuts(int i);
        bool photonPassesMediumIsoCuts(int i);
        bool photonPassesTightIsoCuts(int i);
        bool isTightPhoton(int i);
	double getPhotonScaleCorrectionRunOne7TeV(int run, double eta, double r9, double et);
	double getPhotonScaleCorrectionRunOne8TeV(int run, double eta, double r9, double et);

	//function in HggRazorAuxPhoton.cc
	// R u n 1   C u t   B a s e d   I D
	//----------------------------------
	bool isGoodPhotonRun1( int i, bool _iso, bool _debug );
	bool photonPassIsoRun1( int i, bool _debug );
	void getPhotonEffAreaRun1( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );
	bool passEleVetoRun1( int i );
	// R u n 2   C u t   B a s e d   I D
	//----------------------------------
	bool isGoodPhotonRun2( int i, bool _iso, WP wp, bool _debug );
	bool photonPassIsoRun2( int i, WP wp ,bool _debug );
	void getPhotonEffAreaRun2( float eta, double& effAreaChHad, double& effAreaNHad, double& effAreaPho );


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
	double JetEnergySmearingFactor( double jetPt, double jetEta, double NPU, 
  					SimpleJetResolution *JetResolutionCalculator, TRandom3 *random);
	
        //functions in RazorAuxMisc.cc
	double deltaPhi(double phi1, double phi2);
	double deltaR(double eta1, double phi1, double eta2, double phi2);
        TLorentzVector makeTLorentzVector(double pt, double eta, double phi, double energy);
	TLorentzVector makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass);
	vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
	std::vector< std::vector<int> > getHemispheresV2( std::vector<TLorentzVector> jets);
	double computeMR(TLorentzVector hem1, TLorentzVector hem2);
        double computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector met);
	double GetMT( TLorentzVector visible, TVector3 met );
	double GetMTEnergy( TLorentzVector visible, TVector3 met );
	double GetMT( TLorentzVector visible, TLorentzVector met );
	double GetMTEnergy( TLorentzVector visible, TLorentzVector met );
	double GetDphi( TLorentzVector visible, TVector3 met );
	double GetDphi( TLorentzVector visible, TLorentzVector met );

        bool passesHadronicRazorBaseline(double MR, double Rsq);
        bool passesLeptonicRazorBaseline(double MR, double Rsq);
        int SubtractParticleFromCollection(TLorentzVector ToSubtract, vector<TLorentzVector>& Collection, float deltaRMatch=0.4);
	
	//functions in src/RazorAuxGenLevel.cc
	bool matchesGenMuon(double eta, double phi);
	bool matchesGenElectron(double eta, double phi);
	bool isGenTau(int index);
	bool isGenLeptonicTau(int index);
	int findClosestGenElectron(double eta, double phi);
	int findClosestGenMuon(double eta, double phi);
	int findClosestGenTau(double eta, double phi);
	int findClosestRecoTau(double eta, double phi);
	int GetTauMatchedID(double eta, double phi);
	int findClosestParton(float eta, float phi);

	//Added to src/RazorAuxGenLevel.cc
	int findClosestGenJet(double eta, double phi);
	
        //enums
	// OLD Categories without 6jet category
        /* enum RazorBox { //boxes for razor inclusive analysis */
	/*   MuEle = 0,  */
	/*   MuMu = 1, */
	/*   EleEle = 2, */
	/*   MuMultiJet = 3, */
	/*   MuJet = 4, */
	/*   EleMultiJet = 5, */
	/*   EleJet = 6, */
	/*   LooseLeptonMultiJet = 7, */
	/*   MultiJet = 8, */
	/*   LooseLeptonDiJet = 9, */
	/*   DiJet = 10, */
	/*   TwoBJet = 10, */
	/*   OneBJet = 11, */
	/*   ZeroBJet = 12, */
	/*   NONE = 999 */
        /* }; */
        enum RazorBox { //boxes for razor inclusive analysis
	  MuEle = 0, 
	  MuMu = 1,
	  EleEle = 2,
	  MuSixJet = 3,
	  MuFourJet = 4,
	  MuJet = 5,
	  EleSixJet = 6,
	  EleFourJet = 7,
	  EleJet = 8,
	  LooseLeptonSixJet = 9,
	  LooseLeptonFourJet = 10,
	  SixJet = 11,
	  FourJet = 12,
	  LooseLeptonDiJet = 13,
	  DiJet = 14,	  
	  TwoBJet = 15,
	  OneBJet = 16,
	  ZeroBJet = 17,
	  MuMultiJet = 18,
	  EleMultiJet = 19,
	  LooseLeptonMultiJet = 20,
	  MultiJet = 21,
	  NONE = 999
        };
};

#endif
