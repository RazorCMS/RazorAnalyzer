//Class for analyzing ntuples produced by the RazorTuplizer framework
//
//Author: Caltech Razor team

#ifndef RazorAnalyzer_h
#define RazorAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <RazorEvents.h> //This is a MakeClass of the RazorEvents tree in the ntuple to be analyzed

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

};

#endif
