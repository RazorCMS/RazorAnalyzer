#ifndef ControlSampleEvents_H
#define ControlSampleEvents_H

#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include <cmath>
#include "assert.h"
#include <Rtypes.h>
#include "TLorentzVector.h"

  class ControlSampleEvents {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== Process IDs  ====
      //*******************************************
      enum BkgProcessId { kData = 0,
			  kQCD = 1,
			  kWJets = 2,
			  kZJets = 3,
			  kTTJets = 4,
			  kSingleT = 5,
			  kVV = 6,
			  kHiggs = 7,
			  kSUSY = 99,
			  kUnknown = 999
      };

      /// variables
      Float_t                 weight;
      UInt_t                  run;
      UInt_t                  lumi;
      UInt_t                  event;
      UInt_t                  processID;
      UInt_t                  NPU_0;
      UInt_t                  NPU_Minus1;
      UInt_t                  NPU_Plus1;
      TLorentzVector          genlep1;
      TLorentzVector          genlep2;
      Int_t                   genlep1Type;
      Int_t                   genlep2Type;
      TLorentzVector          lep1;
      TLorentzVector          lep2;
      Int_t                   lep1Type;
      Int_t                   lep2Type;
      Int_t                   lep1MatchedGenLepIndex;
      Int_t                   lep2MatchedGenLepIndex;
      Bool_t                  lep1PassVeto;
      Bool_t                  lep1PassLoose;
      Bool_t                  lep1PassTight;
      Bool_t                  lep2PassVeto;
      Bool_t                  lep2PassLoose;
      Bool_t                  lep2PassTight;
      TLorentzVector          bjet1;
      TLorentzVector          bjet2;
      Bool_t                  bjet1PassLoose;
      Bool_t                  bjet1PassMedium;
      Bool_t                  bjet1PassTight;
      Bool_t                  bjet2PassLoose;
      Bool_t                  bjet2PassMedium;
      Bool_t                  bjet2PassTight;
      Float_t                 MR;
      Float_t                 Rsq;
      Float_t                 MET;
      UInt_t                  NJets40;
      UInt_t                  NBJetsLoose;
      UInt_t                  NBJetsMedium;
      UInt_t                  NBJetsTight;
      Float_t                 HT;
      Float_t                 lep1MT;
      

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      ControlSampleEvents()  {
	genlep1Ptr  = &genlep1;
	genlep2Ptr  = &genlep2;
	lep1Ptr     = &lep1;
	lep2Ptr     = &lep2;
	bjet1Ptr    = &bjet1;
	bjet2Ptr    = &bjet2;       
      };

      /// default destructor
      ~ControlSampleEvents(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
	weight               = 0.0;
	run                  = 0.0;
	lumi                 = 0.0;
	event                = 0.0;
	processID            = ControlSampleEvents::kUnknown;
	NPU_0                = 0.0;
	NPU_Minus1           = 0.0;
	NPU_Plus1            = 0.0;
	genlep1              = TLorentzVector();
	genlep2              = TLorentzVector();
	genlep1Type          = 0.0;
	genlep2Type          = 0.0;
	lep1                 = TLorentzVector();
	lep2                 = TLorentzVector();
	lep1Type             = 0.0;
	lep2Type             = 0.0;
	lep1MatchedGenLepIndex = -1;
	lep2MatchedGenLepIndex = -1;
	lep1PassVeto         = 0.0;
	lep1PassLoose        = 0.0;
	lep1PassTight        = 0.0;
	lep2PassVeto         = 0.0;
	lep2PassLoose        = 0.0;
	lep2PassTight        = 0.0;
	bjet1                = TLorentzVector();
	bjet2                = TLorentzVector();
	bjet1PassLoose       = 0.0;
	bjet1PassMedium      = 0.0;
	bjet1PassTight       = 0.0;
	bjet2PassLoose       = 0.0;
	bjet2PassMedium      = 0.0;
	bjet2PassTight       = 0.0;
	MR                   = 0.0;
	Rsq                  = 0.0;
	MET                  = 0.0;
	NJets40              = 0.0;
	NBJetsLoose          = 0.0;
	NBJetsMedium         = 0.0;
	NBJetsTight          = 0.0;
	HT                   = 0.0;      
	lep1MT               = 0.0;      
      }
    
      /// load a ControlSampleEvents
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("ControlSampleEvent"));
	InitTree();
        assert(tree_);
      }
    
      /// create a ControlSampleEvents
      void CreateTree(){
        tree_ = new TTree("ControlSampleEvent","ControlSampleEvent");
        f_ = 0;

        //book the branches
	tree_->Branch("weight",&weight,"weight/F");
	tree_->Branch("run",&run,"run/i");
	tree_->Branch("lumi",&lumi,"lumi/i");
	tree_->Branch("event",&event,"event/i");
	tree_->Branch("processID",&processID,"processID/i");
	tree_->Branch("NPU_0",&NPU_0,"NPU_0/i");
	tree_->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");
	tree_->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");
	tree_->Branch("genlep1Type",&genlep1Type,"genlep1Type/I");
	tree_->Branch("genlep2Type",&genlep2Type,"genlep2Type/I");
	tree_->Branch("lep1Type",&lep1Type,"lep1Type/I");
	tree_->Branch("lep2Type",&lep2Type,"lep2Type/I");
	tree_->Branch("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex,"lep1MatchedGenLepIndex/I");
	tree_->Branch("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex,"lep2MatchedGenLepIndex/I");
	tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
	tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
	tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
	tree_->Branch("lep2PassVeto",&lep2PassVeto,"lep2PassVeto/O");
	tree_->Branch("lep2PassLoose",&lep2PassLoose,"lep2PassLoose/O");
	tree_->Branch("lep2PassTight",&lep2PassTight,"lep2PassTight/O");
	tree_->Branch("bjet1PassLoose",&bjet1PassLoose,"bjet1PassLoose/O");
	tree_->Branch("bjet1PassMedium",&bjet1PassMedium,"bjet1PassMedium/O");
	tree_->Branch("bjet1PassTight",&bjet1PassTight,"bjet1PassTight/O");
	tree_->Branch("bjet2PassLoose",&bjet2PassLoose,"bjet2PassLoose/O");
	tree_->Branch("bjet2PassMedium",&bjet2PassMedium,"bjet2PassMedium/O");
	tree_->Branch("bjet2PassTight",&bjet2PassTight,"bjet2PassTight/O");
	tree_->Branch("MR",&MR,"MR/F");
	tree_->Branch("Rsq",&Rsq,"Rsq/F");
	tree_->Branch("MET",&MET,"MET/F");
	tree_->Branch("NJets40",&NJets40,"NJets40/i");
	tree_->Branch("NBJetsLoose",&NBJetsLoose,"NBJetsLoose/i");
	tree_->Branch("NBJetsMedium",&NBJetsMedium,"NBJetsMedium/i");
	tree_->Branch("NBJetsTight",&NBJetsTight,"NBJetsTight/i");
	tree_->Branch("HT",&HT,"HT/F");
	tree_->Branch("lep1MT",&lep1MT,"lep1MT/F");
	tree_->Branch("genlep1", "TLorentzVector", &genlep1Ptr);
	tree_->Branch("genlep2", "TLorentzVector", &genlep2Ptr);
	tree_->Branch("lep1",    "TLorentzVector", &lep1Ptr);
	tree_->Branch("lep2",    "TLorentzVector", &lep2Ptr);
	tree_->Branch("bjet1",   "TLorentzVector", &bjet1Ptr);
	tree_->Branch("bjet2",   "TLorentzVector", &bjet2Ptr);
      } 

      // initialze a ControlSampleEvents
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

	tree_->SetBranchAddress("weight",&weight);
	tree_->SetBranchAddress("run",&run);
	tree_->SetBranchAddress("lumi",&lumi);
	tree_->SetBranchAddress("event",&event);
	tree_->SetBranchAddress("processID",&processID);
	tree_->SetBranchAddress("NPU_0",&NPU_0);
	tree_->SetBranchAddress("NPU_Minus1",&NPU_Minus1);
	tree_->SetBranchAddress("NPU_Plus1",&NPU_Plus1);
	tree_->SetBranchAddress("genlep1Type",&genlep1Type);
	tree_->SetBranchAddress("genlep2Type",&genlep2Type);
	tree_->SetBranchAddress("lep1Type",&lep1Type);
	tree_->SetBranchAddress("lep2Type",&lep2Type);
	tree_->SetBranchAddress("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex);
	tree_->SetBranchAddress("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex);
	tree_->SetBranchAddress("lep1PassVeto",&lep1PassVeto);
	tree_->SetBranchAddress("lep1PassLoose",&lep1PassLoose);
	tree_->SetBranchAddress("lep1PassTight",&lep1PassTight);
	tree_->SetBranchAddress("lep2PassVeto",&lep2PassVeto);
	tree_->SetBranchAddress("lep2PassLoose",&lep2PassLoose);
	tree_->SetBranchAddress("lep2PassTight",&lep2PassTight);
	tree_->SetBranchAddress("bjet1PassLoose",&bjet1PassLoose);
	tree_->SetBranchAddress("bjet1PassMedium",&bjet1PassMedium);
	tree_->SetBranchAddress("bjet1PassTight",&bjet1PassTight);
	tree_->SetBranchAddress("bjet2PassLoose",&bjet2PassLoose);
	tree_->SetBranchAddress("bjet2PassMedium",&bjet2PassMedium);
	tree_->SetBranchAddress("bjet2PassTight",&bjet2PassTight);
	tree_->SetBranchAddress("MR",&MR);
	tree_->SetBranchAddress("Rsq",&Rsq);
	tree_->SetBranchAddress("MET",&MET);
	tree_->SetBranchAddress("NJets40",&NJets40);
	tree_->SetBranchAddress("NBJetsLoose",&NBJetsLoose);
	tree_->SetBranchAddress("NBJetsMedium",&NBJetsMedium);
	tree_->SetBranchAddress("NBJetsTight",&NBJetsTight);
	tree_->SetBranchAddress("HT",&HT);
	tree_->SetBranchAddress("lep1MT",&lep1MT);	
	tree_->SetBranchAddress("genlep1",&genlep1Ptr);
	tree_->SetBranchAddress("genlep2",&genlep2Ptr);
	tree_->SetBranchAddress("lep1",&lep1Ptr);
	tree_->SetBranchAddress("lep2",&lep2Ptr);
	tree_->SetBranchAddress("bjet1",&bjet1Ptr);
	tree_->SetBranchAddress("bjet2",&bjet2Ptr);
	
        gErrorIgnoreLevel = currentState;
      }

    private:
      TLorentzVector* genlep1Ptr;
      TLorentzVector* genlep2Ptr;
      TLorentzVector* lep1Ptr;
      TLorentzVector* lep2Ptr;
      TLorentzVector* bjet1Ptr;
      TLorentzVector* bjet2Ptr;
      
  }; 


#endif

