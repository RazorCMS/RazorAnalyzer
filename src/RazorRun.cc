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
    if(argc < 3){
        cerr << "Usage: RazorRun <input list> <analysis type> <output filename (optional)> [option number] [optional label]" << endl;
        cerr << "Analyses available: " << endl 
            << "razor                --   inclusive razor analysis" << endl 
            << "hggrazor             --   higgs->diphoton razor analysis" << endl
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

    string outputFileName = "";
    if (argc >= 4)  outputFileName = argv[3];

    int option = -1;
    if (argc >= 5) {
      option = atoi(argv[4]);
    }
    
    string label = "";
    if (argc >= 6) {
      label = argv[5];
    }
    
    //build the TChain
    //TChain *theChain = new TChain("ntuples/RazorEvents");
    TChain* theChain = new TChain("RazorEvents"); 
    string curFileName;
    ifstream inputFile(inputFileName.c_str());
    int NFilesLoaded = 0;
    if(!inputFile){
        cerr << "Error: input file not found!" << endl; 
        return -1;
    }
    while(getline(inputFile, curFileName)){
        theChain->Add(curFileName.c_str());
	if (analysisType != "MakeMCPileupDistribution") {
	  std::cout << "chaining " << curFileName << std::endl;
	}
	NFilesLoaded++;
    }
    std::cout << "Loaded Total of " << NFilesLoaded << " files\n";

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
	analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
	analyzer.EnableMC();
	analyzer.EnableGenParticles();
        analyzer.RazorInclusive(outputFileName, true); //change the bool to true if you want all analysis boxes combined in one tree
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
    else if(analysisType == "photonntupler"){
        cout << "Running photon ntupler..." << endl;
        analyzer.EnableEventInfo();
        analyzer.EnablePhotons();
        analyzer.EnableGenParticles();
        analyzer.PhotonNtupler();
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
      analyzer.RazorControlRegions(outputFileName, option, false);
    }
    else if(analysisType == "RunOneRazorControlRegions"){
      cout << "Executing RunOneRazorControlRegions analysis..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();      
      analyzer.RazorControlRegions(outputFileName, option, true);
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
      if(option == 1){ 
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
      if(option == 1){ 
          analyzer.RazorPhotonStudy(outputFileName, true, false, true); //run with data -- don't filter events (Run 1)
      }
      else{
          analyzer.RazorPhotonStudy(outputFileName, false, false, true); //run with MC -- don't filter events (Run 1)
      }
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
