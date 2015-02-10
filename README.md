RazorAnalyzer
=============

Class for analyzing the 2015 razor ntuples

Setup
-------------
From t3-higgs:

    git clone https://github.com/RazorCMS/RazorAnalyzer.git
    cd RazorAnalyzer
    make
  
Defining a new analysis
-------------
1) Copy analyses/DummyAnalysis.cc and modify it to define your analyzer's behavior.  Be sure to change the name of the DummyAnalysis() function.

2) Add your analysis under "LIST OF ANALYSES" in include/RazorAnalyzer.h.

3) In src/RazorRun.cc, under "EXECUTE YOUR ANALYSIS" add the option to execute your analysis code.

Running
------------
After compiling, 

    ./RazorRun <list of input files> <analysis type>
  
Example: to execute a dummy analysis that does nothing,

    ./RazorRun lists/TTJets_List_Test.txt dummy

Normalizing the processed ntuples
------------
The NormalizeNtuple macro opens a specified set of files and adds a 'weight' branch to each TTree in each file.  The value of 'weight' is the same for all events in a tree and is equal to CrossSection/NEvents, where NEvents is the total number of events processed for the given dataset.  The cross sections can be found in the file data/xSections.dat.  To run NormalizeNtuple:

    ./NormalizeNtuple <input file>

See lists/filestonormalize/testTTJets.txt for an example input file to be used with NormalizeNtuple.

Running over many samples at once
-----------
The script 'runEverything' takes as input the name of a directory containing lists of files to process.  Ex:

    ./runEverything lists/27Oct2014/

The script will (locally) run the RazorInclusive analysis on each sample listed in the given directory, run NormalizeNtuple to add event weights, and add all of the output files together.

Fitting samples
-----------
Make directories

	mkdir Backgrounds; mkdir Signals; mkdir Datasets

Copy background trees such as
RazorAnalysis\_TTJets\_25ns\_1pb\_weighted.root to Backgrounds/ and signal
trees such as
RazorAnalysis\_SMS-T1bbbb\_2J_mGl-1500\_mLSP-100\_25ns\_1pb\_weighted.root
to Signals/. Note scripts assume no QCD and scale up all backgrounds
by hard-coded factors per box.

	mkdir Backgrounds; mkdir Signals; mkdir Datasets
	python python/DustinTuple2RooDataSet.py -b MultiJet -c config/run2.config -w -l 4000 -d Datasets Backgrounds/RazorAnalysis_*_25ns_1pb_weighted.root 
	python python/RooDataSet2UnweightedDataSet.py -b MultiJet -c config/run2.config -d Datasets Datasets/RazorAnalysis_SMCocktail_weighted_lumi-4.0_MultiJet.root
	python python/Fit.py -b MultiJet -c config/run2.config -d FitResults Datasets/RazorAnalysis_SMCocktail_unweighted_lumi-4.0_MultiJet.root
	
The scripts will (locally) produce the weighted and "unweighted" datasets, run the fit, and make
the 1D fit projections in MR and Rsq.


