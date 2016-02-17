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

std::string ParseCommandLine( int argc, char* argv[], std::string opt )
{
  for (int i = 1; i < argc; i++ )
    {
      std::string tmp( argv[i] );
      if ( tmp.find( opt ) != std::string::npos )
        {
          if ( tmp.find( "=" )  != std::string::npos ) return tmp.substr( tmp.find_last_of("=") + 1 );
	  if ( tmp.find( "--" ) != std::string::npos ) return "yes";
	}
    }
  
  return "";
};

void usage()
{
  std::cerr << "Usage: RazorRun  <input list>  <analysis type>  [options]\n[options]:\n"
	    << "-d  --isData\n"
	    << "-f  --outputFile=<output filename> (optional)\n"
	    << "-n  --optionNumber=<option number> (optional)\n"
	    << "-l  --optionLabel=<option Label> (optional)\n" 
	    << "-h  --help"
	    << std::endl;
  
  std::cerr << "Analyses available:\n"
	    << "razor                --   inclusive razor analysis\n"
	    << "razorfull            --   larger razor ntuple used for setting limits w/systematics\n"
	    << "hggrazor             --   higgs->diphoton razor analysis\n"
	    << "hggrazorExo15004     --   higgs->diphoton razor analysis using Exo15004\n"
	    << "matchedrazor         --   inclusive razor analysis using only jets matched to genjets\n"
	    << "razorVetoLeptonStudy --   study lepton veto\n"
	    << "razorPhotonStudy     --   select events for Z->invisible control sample\n"
	    << "razorPhotonStudyEff  --   same tree content as razorPhotonStudy, but save all events\n"
	    << "electronNtupler      --   study electron variables\n"
	    << "muonNtupler          --   study muon variables\n"
	    << "jetNtupler           --   study jet variables\n"
	    << "photonntupler        --   study photon variables\n"
	    << "dummy                --   do nothing useful\n"
	    << "RazorDM              --   run MultiJet razor dark matter analysis\n"
            << "RazorPhotonDM        --   run photon+jets razor dark matter analysis\n"
        << "RazorAlphaT          --   run MultiJet alphaT analysis\n"
	    << "RazorTagAndProbe     --   run tag and probe using Z peak" 
	    << std::endl;
};

int main(int argc, char* argv[]){

  //get input files and analysis type from command line
  if ( ParseCommandLine( argc, argv, "--help" ) != ""  || ParseCommandLine( argc, argv, "-h" ) != ""  || argc < 3 )
    {
      usage();
      return -1;
    }
  
  //----------------------------------------
  //Getting <input list> and <analysis type>
  //----------------------------------------
  string inputFileName(argv[1]);
  string analysisType(argv[2]);
  
  //--------------------------------
  //G e t t i n g   d a t a  f l a g 
  //--------------------------------
  std::string _isData = ParseCommandLine( argc, argv, "--isData" );
  std::string _d = ParseCommandLine( argc, argv, "-d" );
  bool isData = false;
  if ( _isData == "yes" || _d == "yes" ) isData = true;

  //---------------------------------------------
  //G e t t i n g   o u t p u t F i l e   N a m e  
  //---------------------------------------------
  std::string _outFile = ParseCommandLine( argc, argv, "--outputFile=" );
  std::string _f = ParseCommandLine( argc, argv, "-f=" );
  string outputFileName = "";
  if ( _outFile != "" )
    {
      outputFileName = _outFile;
    }
  else if ( _f != "" )
    {
      outputFileName = _f;
    }
  else
    {
      std::cerr << "[WARNING]: output ROOT file not provided, using default output" << std::endl;
    }
  
  //-----------------------------------------
  //G e t t i n g   o p t i o n   n u m b e r
  //-----------------------------------------
  int option = -1;
  std::string _optionNumber = ParseCommandLine( argc, argv, "--optionNumber=" );
  std::string _n = ParseCommandLine( argc, argv, "-n=" );
  if ( _optionNumber != "" )
    {
      option = atoi( _optionNumber.c_str() );
    }
  else if ( _n != "" )
    {
      option = atoi( _n.c_str() );
    } 
  else
    {
      std::cerr << "[WARNING]: option number not provided, using default option number" << std::endl;
    }
  
  string label = "";
  std::string _optionLabel = ParseCommandLine( argc, argv, "--optionLabel=" );
  std::string _l = ParseCommandLine( argc, argv, "-l=" );
  if ( _optionLabel != "" ) 
    {
      label = _optionLabel;
    }
  else if ( _l != "" )
    {
      label = _l;
    }
  else
    {
      std::cerr << "[WARNING]: optional label not provided, using default optional label" << std::endl;
    }
  
  std::cout << "[INFO]: <input list> --> " << inputFileName << std::endl;
  std::cout << "[INFO]: <analysis type> --> " << analysisType << std::endl;
  std::cout << "[INFO]: isData --> " << isData << std::endl;
  std::cout << "[INFO]: outputFileName --> " << outputFileName << std::endl;
  std::cout << "[INFO]: option --> " << option << std::endl;
  std::cout << "[INFO]: optionalLabel --> " << label << std::endl;
    
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
    else if(analysisType == "razor" || analysisType == "RazorInclusive"){
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
        bool isFastsimSMS = (option == 1); //specify option = 1 to split fastsim signal samples by mass point
        analyzer.RazorInclusive(outputFileName, true, isData, isFastsimSMS); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "razorfull" || analysisType == "fullrazor" || analysisType == "FullRazorInclusive"){
        cout << "Executing full razor inclusive analysis..." << endl;
        analyzer.EnableAll();
        bool isFastsimSMS = (option == 1); //specify option = 1 to split fastsim signal samples by mass point
        analyzer.FullRazorInclusive(outputFileName, isData, isFastsimSMS); 
    }
    else if(analysisType == "hggrazor" || analysisType == "HggRazor"){
        cout << "Executing higgs->diphoton razor analysis..." << endl;
        analyzer.EnableEventInfo();
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
        analyzer.EnablePhotons();
	analyzer.EnableMC();
	analyzer.EnableGenParticles();
	analyzer.EnablePileup();      
        analyzer.HggRazor(outputFileName, true, option, isData); //change the bool to true if you want all analysis boxes combined in one tree
    }
    else if ( analysisType == "hggrazorgen" )
      {
	cout << "Executing higgs->diphoton razor analysis GEN LEVEL..." << endl;
        analyzer.EnableEventInfo();
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
        analyzer.EnableTaus();
        analyzer.EnablePhotons();
        analyzer.EnableMC();
        analyzer.EnableGenParticles();
        analyzer.EnablePileup();
	analyzer.HggRazorGenLevel(outputFileName, true, option, isData);
      }
    else if(analysisType == "hggrazorExo15004" || analysisType == "hggrazorexo15004"){
      cout << "Executing higgs->diphoton razor analysis (EXO15004)..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnablePhotons();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();
      analyzer.HggRazorExo15004(outputFileName, true, option, isData);
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
    else if(analysisType == "razorVetoLeptonStudy" || analysisType == "RazorVetoLeptonStudy"){
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
    else if(analysisType == "razorZAnalysis" || analysisType == "RazorZAnalysis"){
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
    else if(analysisType == "electronNtupler" || analysisType == "ElectronNtupler"){
      cout << "Executing electron ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableElectrons();
      analyzer.EnableGenParticles();
      analyzer.ElectronNtupler(outputFileName, option);
    }
    else if(analysisType == "muonNtupler" || analysisType == "MuonNtupler"){
      cout << "Executing muon ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableMuons();
      analyzer.EnableGenParticles();
      analyzer.MuonNtupler(outputFileName, option);
    }
    else if(analysisType == "tauNtupler" || analysisType == "TauNtupler"){
      cout << "Executing tau ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableTaus();
      analyzer.EnableGenParticles();
      analyzer.TauNtupler(outputFileName, option);
    }
    else if(analysisType == "jetNtupler" || analysisType == "JetNtupler"){
      cout << "Executing jet ntupler..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableGenParticles();
      analyzer.EnableMC();
      analyzer.JetNtupler(outputFileName, option);
    }
    else if(analysisType == "photonNtupler" || analysisType == "PhotonNtupler"){
        cout << "Running photon ntupler..." << endl;
        analyzer.EnableEventInfo();
        analyzer.EnablePhotons();
        analyzer.EnableGenParticles();
        analyzer.PhotonNtupler(outputFileName, option);
    }
    else if(analysisType == "met" || analysisType == "RazorMetAna"){ // met analyzer to plot some histograms
        cout << "Executing razor MET analysis..." << endl;
        analyzer.EnableJets();
        analyzer.EnableMet();
        analyzer.EnableElectrons();
        analyzer.EnableMuons();
	    analyzer.EnableEventInfo();
        analyzer.RazorMetAna(outputFileName); //change to true if you want all analysis boxes combined in one tree
    }
    else if(analysisType == "RazorDM" || analysisType == "razorDM"){
      cout << "Executing RazorDM analysis..." << endl;
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableEventInfo();
      analyzer.EnableMC();
      analyzer.RazorDM(outputFileName, true, isData);
    }else if(analysisType == "RazorPhotonDM" || analysisType == "razorPhotonDM"){
      cout << "Executing RazorPhotonDM analysis..." << endl;
      analyzer.EnablePhotons();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableEventInfo();
      analyzer.EnableMC();
      analyzer.RazorPhotonDM(outputFileName, true, isData);
    }else if(analysisType == "RazorAlphaT" || analysisType == "razorAlphaT"){
      cout << "Executing RazorAlphaT analysis..." << endl;
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnableEventInfo();
      analyzer.EnableMC();
      analyzer.RazorAlphaT(outputFileName, true, isData);
    }
    else if(analysisType == "RazorControlRegions"){
      cout << "Executing RazorControlRegions analysis..." << endl;
      analyzer.EnableEventInfo();
      analyzer.EnableJets();
      analyzer.EnableMet();
      analyzer.EnableElectrons();
      analyzer.EnableMuons();
      analyzer.EnableTaus();
      analyzer.EnablePhotons();
      analyzer.EnableMC();
      analyzer.EnableGenParticles();
      analyzer.EnablePileup();      
      analyzer.RazorControlRegions(outputFileName, option, isData);
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
    else if(analysisType == "razorPhotonStudy" || analysisType == "RazorPhotonStudy"){
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
    else if(analysisType == "hbbrazor" || analysisType == "HbbRazor"){
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
    else if(analysisType == "hzzRazor" || analysisType == "HZZRazor"){
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
        analyzer.HZZRazor(outputFileName, isData); //change the bool to true if you want all analysis boxes combined in one tree
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
      analyzer.EnableIsoPFCandidates();
      bool isRunOne = false; 
      if (analysisType == "RazorRunOneQCDStudy") { isRunOne = true; }
      analyzer.RazorQCDStudy(outputFileName, isData, isRunOne); 
    }
    else if(analysisType == "RazorTagAndProbe"){
      cout << "Executing RazorTagAndProbe analysis..." << endl;
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
      analyzer.RazorTagAndProbe(outputFileName, option, isData);
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
