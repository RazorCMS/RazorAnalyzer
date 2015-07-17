//C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>

#include <RazorAnalyzer.h>

using namespace std;

int main(int argc, char* argv[]){

    //get input files and analysis type from command line
    if(argc < 4){
        cerr << "Usage: RazorRun <input list> <analysis type> <isData> <output filename (optional)> [option number] [optional label]" << endl;
        cerr << "Analyses available: " << endl 
	     << "razor                --   inclusive razor analysis" << endl 
	     << "hggrazor             --   run1 higgs->diphoton razor analysis" << endl
	     << "run2hggrazor         --   run2 higgs->diphoton razor analysis" << endl
	     << "matchedrazor         --   inclusive razor analysis using only jets matched to genjets" << endl 
	     << "razorVetoLeptonStudy --   study lepton veto" << endl
	     << "razorPhotonStudy     --   select events for Z->invisible control sample" << endl
	     << "razorPhotonStudyEff  --   same tree content as razorPhotonStudy, but save all events" << endl
	     << "electronNtupler      --   study electron variables" << endl
	     << "muonNtupler          --   study muon variables" << endl
	     << "jetNtupler           --   study jet variables" << endl
	     << "photonntupler        --   study photon variables" << endl
	     << "dummy                --   do nothing useful" << endl;
        return -1;
    }
    string inputFileName(argv[1]);
    string analysisType(argv[2]);

    string isDataOption = "";
    isDataOption = argv[3];
    bool isData = false;
    if (isDataOption == "y" || isDataOption == "1" || isDataOption == "true") isData = true;
        
    string outputFileName = "";
    if (argc >= 5)  outputFileName = argv[4];

    int option = -1;
    if (argc >= 6) {
      option = atoi(argv[5]);
    }
    
    string label = "";
    if (argc >= 7) {
      label = argv[6];
    }
    
    //build the TChain
    //tree name is set give the structure in the first root file, see while loop below
    TChain* theChain = new TChain();
    string curFileName;
    ifstream inputFile(inputFileName.c_str());
    int NFilesLoaded = 0;
    if ( !inputFile ){
      cerr << "Error: input file not found!" << endl;
      return -1;
    }
    
    while ( getline(inputFile, curFileName) )
      {
	if ( NFilesLoaded == 0 )
	  {
	    /*
	      checks root file structure and add first file
	    */
	    TFile* f_0 = TFile::Open( curFileName.c_str() );
	    if( f_0->GetDirectory("ntuples") )
	      {
		theChain->SetName("ntuples/RazorEvents");
		std::cout << "[INFO]: default configuration for tchain" << std::endl;
	      }
	    else
	      {
		theChain->SetName("RazorEvents");
		std::cout << "[INFO]: alternative configuration for tchain"<< std::endl;
	      }
	    theChain->Add( curFileName.c_str() );
	    delete f_0;
	  }
	else
	  {
	    //Addind remaining files after file structure is decided
	    theChain->Add( curFileName.c_str() );
	  }
	
	if ( analysisType != "MakeMCPileupDistribution" ) {
	  std::cout << "chaining " << curFileName << std::endl;
        }
        NFilesLoaded++;
      }     
    
    std::cout << "Loaded Total of " << NFilesLoaded << " files\n";
    if ( theChain == NULL ) return -1;
    
    RazorAnalyzer analyzer(theChain);
    
    //------ EXECUTE YOUR ANALYSIS ------//
    if(analysisType == "dummy"){
        cout << "Executing dummy analysis..." << endl;
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.DummyAnalysis();
    }
    else if(analysisType == "razor"){
        cout << "Executing razor inclusive analysis..." << endl;
	analyzer.EnableEventInfo();
	analyzer.EnablePileup();
	analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
	analyzer.EnableMC();
	analyzer.EnableGenParticles();
        analyzer.RazorInclusive(outputFileName, true, isData, false); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "RunOneRazor"){
        cout << "Executing razor inclusive analysis..." << endl;
	analyzer.EnableEventInfo();
	analyzer.EnablePileup();
	analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
	analyzer.EnableMC();
	analyzer.EnableGenParticles();
        analyzer.RazorInclusive(outputFileName, true, isData, true); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "hggrazor"){
        cout << "Executing higgs->diphoton razor analysis..." << endl;
        analyzer.EnableEventInfo();
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
        analyzer.EnablePhotons();
        analyzer.HggRazor(outputFileName, false); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "run2hggrazor"){
      cout << "Executing run2 higgs->diphoton razor analysis..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnablePhotons();
      analyzer.HggRazorRun2(outputFileName, false); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "matchedrazor"){
      cout << "Executing genjet-matched razor inclusive analysis..." << endl;
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
        analyzer.EnableMC();
        analyzer.MatchedRazorInclusive(outputFileName, false); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "razorVetoLeptonStudy"){
      cout << "Executing razorVetoLeptonStudy..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();
      analyzer.EnableIsoPFCandidates();
      analyzer.EnablePhotons();
      if (option == 1) {
	analyzer.RazorVetoLeptonStudy(outputFileName, true);
      } else {
	analyzer.RazorVetoLeptonStudy(outputFileName, false);
      }
    }
    else if(analysisType == "razorZAnalysis"){
      cout << "Executing razorZAnalysis..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();
      analyzer.EnableIsoPFCandidates();
      analyzer.RazorZAnalysis(outputFileName, true);
    }
    else if(analysisType == "electronNtupler"){
      cout << "Executing electron ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableElectrons();
      analyzer.EnableGenParticles();
      analyzer.ElectronNtupler(outputFileName, option);
    }
    else if(analysisType == "muonNtupler"){
      cout << "Executing muon ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableMuons();
      analyzer.EnableGenParticles();
      analyzer.MuonNtupler(outputFileName, option);
    }
    else if(analysisType == "tauNtupler"){
      cout << "Executing tau ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableTaus();
      analyzer.EnableGenParticles();
      analyzer.TauNtupler(outputFileName, option);
    }
    else if(analysisType == "jetNtupler"){
      cout << "Executing jet ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableGenParticles();
      analyzer.EnableMC();
      analyzer.JetNtupler(outputFileName, option);
    }
    else if(analysisType == "photonNtupler"){
        cout << "Running photon ntupler..." << endl;
        analyzer.EnableEventInfo();
        analyzer.EnablePhotons();
        analyzer.EnableGenParticles();
        analyzer.PhotonNtupler(outputFileName, option);
    }
    else if(analysisType == "met"){ // met analyzer to plot some histograms
        cout << "Executing razor MET analysis..." << endl;
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
	analyzer.EnableEventInfo();
        analyzer.RazorMetAna(outputFileName); //change to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "RazorDM"){
      cout << "Executing RazorDM analysis..." << endl;
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableEventInfo();
      analyzer.EnableMC();

      analyzer.RazorDM(outputFileName);
    }
    else if(analysisType == "RazorControlRegions"){
      cout << "Executing RazorControlRegions analysis..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();      
      analyzer.RazorControlRegions(outputFileName, option, isData, false);
    }
    else if(analysisType == "RunOneRazorControlRegions"){
      cout << "Executing RunOneRazorControlRegions analysis..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnablePhotons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();      
      analyzer.RazorControlRegions(outputFileName, option, isData, true);
    }
    else if(analysisType == "VetoLeptonEfficiencyControlRegion"){
      cout << "Executing VetoLeptonEfficiencyDileptonControlRegion analysis..." << endl;
      analyzer.EnableEventInfo(); 
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();      
      analyzer.VetoLeptonEfficiencyControlRegion(outputFileName, option);
    }
    else if(analysisType == "razorPhotonStudy"){
      cout << "Executing razorPhotonStudy..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnablePhotons();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();
      if(isData){ 
          analyzer.RazorPhotonStudy(outputFileName, true, true, true); //run with data (Run 1)
      }
      else{
          analyzer.RazorPhotonStudy(outputFileName, false, true, true); //run with MC (Run 1)
      }
    }    
    else if(analysisType == "razorPhotonStudyEff"){
      cout << "Executing razorPhotonStudy, keeping all events for efficiency calculation..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnablePhotons();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();
      if(isData){ 
          analyzer.RazorPhotonStudy(outputFileName, true, false, true); //run with data -- don't filter events (Run 1)
      }
      else{
          analyzer.RazorPhotonStudy(outputFileName, false, false, true); //run with MC -- don't filter events (Run 1)
      }
    } 
    else if(analysisType == "hbbrazor"){
        cout << "Executing razor inclusive analysis..." << endl;
	analyzer.EnableEventInfo();
	analyzer.EnablePileup();
	analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
	analyzer.EnableMC();
	analyzer.EnableGenParticles();
        analyzer.HbbRazor(outputFileName, true, isData, false); 
    }
    else if(analysisType == "hzzRunOneRazor"){
        cout << "Executing razor inclusive analysis..." << endl;
	analyzer.EnableEventInfo();
	analyzer.EnablePileup();
	analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
	analyzer.EnableMC();
	analyzer.EnableGenParticles();
        analyzer.HZZRazor(outputFileName, isData, true); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "hzzRazor"){
        cout << "Executing razor inclusive analysis..." << endl;
	analyzer.EnableEventInfo();
	analyzer.EnablePileup();
	analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
	analyzer.EnableMC();
	analyzer.EnableGenParticles();
        analyzer.HZZRazor(outputFileName, isData, false); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "RazorQCDStudy" || analysisType == "RazorRunOneQCDStudy" ){
      cout << "Executing razor QCD analysis..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnablePileup();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      bool isRunOne = false; 
      if (analysisType == "RazorRunOneQCDStudy") isRunOne = true;
      analyzer.RazorQCDStudy(outputFileName, isRunOne, false); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "MakeMCPileupDistribution"){
      cout << "Executing MakeMCPileupDistribution..." << endl;
      analyzer.EnablePileup();     
      analyzer.MakeMCPileupDistribution(outputFileName,label);
    }    



    else { //analysis not found
      cerr << "Error: the given analysis type is not defined in RazorRun.cc!" << endl;
    }

    cout << "Process completed!" << endl;
    cerr << "------------------------------" << endl; //present so that an empty .err file corresponds to a failed job
    
    return 0;
}
