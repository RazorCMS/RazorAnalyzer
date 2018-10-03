#ifndef DiMuonTree_H
#define DiMuonTree_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"

  class DiMuonTree {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== DiMuonTriggerBits  ====
      //*******************************************
      enum DiMuonTriggerBits { kDiMuonTrigger_DiMuon                                    = 0x000001
      };

      /// variables
      //UInt_t                  fi_evt;
      Float_t                 f_weight;
      UInt_t                  f_run_number;
      UInt_t                  f_lumi_section_number;
      UInt_t                  f_event_number;
      //DiMuon
      Float_t                 f_mumu_mass;
      Float_t                 f_mumu_pt;
      Float_t                 f_mumu_eta;
      Float_t                 f_mumu_phi;
      //Muons
      Float_t                 f_mu_pt[2];
      Float_t                 f_mu_eta[2];
      Float_t                 f_mu_phi[2];
      //
      Float_t                 f_met_pt;
      Float_t                 f_met_phi;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;

      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor
      DiMuonTree()  {};
      /// default destructor
      ~DiMuonTree(){
        if (f_) f_->Close();
      };

      /// initialize varibles and fill list of available variables
      void InitVariables()
      {
        fi_evt = 0;
        ResetVariables();
      }
      void ResetVariables() {
        f_weight			                     = 0.0;
        f_run_number		                     = 0.0;
        f_lumi_section_number	                     = 0.0;
        f_event_number		                     = 0.0;
        //
        f_mumu_mass 	                     = 0.0;
        f_mumu_pt 		                     = 0.0;
        f_mumu_eta 		                     = 0.0;
        f_mumu_phi 		                     = 0.0;
        //
        for ( int i = 0; i < 2; i++ )
        {
          f_mu_pt[i]  = 0.0;
          f_mu_eta[i] = 0.0;
          f_mu_phi[i] = 0.0;
        }
        //
        f_met_pt 		                     = 0.0;
        f_met_phi 			                 = 0.0;
      }

      /// load a loadTree
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("DiMuon"));
        assert(tree_);
      }

      /// create a JetTree
      void CreateTree(){
        tree_ = new TTree("DiMuon","DiMuon");
        f_ = 0;

        //book the branches
        tree_->Branch("weight",&f_weight,"weight/F");
        tree_->Branch("run",&f_run_number,"run/i");
        tree_->Branch("lumi",&f_lumi_section_number,"lumi/i");
        tree_->Branch("event",&f_event_number,"event/i");
        //
        tree_->Branch("mumu_mass",&f_mumu_mass,"mumu_mass/F");
        tree_->Branch("mumu_pt",&f_mumu_pt,"mumu_pt/F");
        tree_->Branch("mumu_eta",&f_mumu_eta,"mumu_eta/F");
        tree_->Branch("mumu_phi",&f_mumu_phi,"mumu_phi/F");
        //
        tree_->Branch("mu_pt",&f_mu_pt,"mu_pt/F");
        tree_->Branch("mu_eta",&f_mu_eta,"mu_eta/F");
        tree_->Branch("mu_phi",&f_mu_phi,"mu_phi/F");
        //
        tree_->Branch("met_pt",&f_met_pt,"met_pt/F");
        tree_->Branch("met_phi",&f_met_phi,"met_phi/F");
      }

      // initialze a JetTree
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

        tree_->SetBranchAddress("weight",&f_weight);
        tree_->SetBranchAddress("run",&f_run_number);
        tree_->SetBranchAddress("lumi",&f_lumi_section_number);
        tree_->SetBranchAddress("event",&f_event_number);
        //
        tree_->Branch("mumu_mass",&f_mumu_mass);
        tree_->Branch("mumu_pt",&f_mumu_pt);
        tree_->Branch("mumu_eta",&f_mumu_eta);
        tree_->Branch("mumu_phi",&f_mumu_phi);
        //
        tree_->Branch("mu_pt",&f_mu_pt);
        tree_->Branch("mu_eta",&f_mu_eta);
        tree_->Branch("mu_phi",&f_mu_phi);
        //
        tree_->Branch("met_pt",&f_met_pt);
        tree_->Branch("met_phi",&f_met_phi);

        gErrorIgnoreLevel = currentState;
      }

  };


#endif
