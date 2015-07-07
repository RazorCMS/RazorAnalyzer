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

Fitting samples and setting limits
-----------

The main controller of the fit and limit setting comes from the
configuration file, which defines the binning used in the binned fit,
the initlial values of the shape parameters, etc.

    config/run2_1-3btag.config

Setup combine from lxplus:

    cmsrel CMSSW_7_1_5
    cd CMSSW_7_1_5/src/
    git clone https://github.com/RazorCMS/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git pull origin razor1dpdf_71X
    scramv1 b clean; scramv1 b
    cd ../..

Proceed with the usual setup of RazorAnalyzer:

    git clone https://github.com/RazorCMS/RazorAnalyzer.git
    cd RazorAnalyzer
    make
  
Make output directories to make workflow easier 

    mkdir Backgrounds; mkdir Signals; mkdir Datasets
    mkdir FitResults; mkdir FitProjections; mkdir cards

Copy all the SM background trees such as RazorInclusive\_TTJets\_1pb\_weighted.root to Backgrounds and signal
trees such as RazorInclusive\_SMS-T1bbbb\_2J\_mGl-1500\_mLSP-100\_1pb\_weighted.root
to Signals.

Note the first script removes the QCD contribution and scales up the
remaining backgrounds if given the option -q. The following commands will (locally)
produce the weighted and "unweighted" datasets, run the fit, and make the 1D fit projections in MR and Rsq. 

    python python/DustinTuple2RooDataSet.py -b MultiJet -c config/run2_1-3btag.config -w -l 3000 -d Datasets Backgrounds/RazorInclusive_*weighted.root
    python python/BinnedFit.py -b MultiJet -c config/run2_1-3btag.config -d FitResults Datasets/RazorInclusive_SMCocktail_weighted_lumi-3.0_1-3btag_MultiJet.root
	
To produce a large number of other 1D fit projections, and the 2D
projections in MR and Rsq,

    python python/PlotFit.py -b MultiJet -c config/run2_1-3btag.config -d FitResults -i FitResults/BinnedFitResults_MultiJet.root
	
Next, to produce the datacards, run combine, and make the fit
projection plots (on lxplus), execute the following:

    python python/DustinTuple2RooDataSet.py -b MultiJet -c config/run2_1-3btag.config -w -l 3000 -d Datasets Signals/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root
    python python/RunCombine.py -b MultiJet -c config/run2_1-3btag.config -d cards --lumi-array 0.2,0.5,1,3,5,7,10 -m T1bbbb --mGluino 1500 --mLSP 100
    python python/PlotCombine.py -b MultiJet -c config/run2_1-3btag.config -d FitProjections -i cards --lumi-array 0.2,0.5,1,3,5,7,10 -m T1bbbb --mGluino 1500 --mLSP 100 
	
which should produce the fit projection plots in the directory
FitProjections.

To make an "unweighted" dataset (not needed for the preceding commands),
execute

    python python/RooDataSet2UnweightedDataSet.py -b MultiJet -c config/run2_1-3btag.config -d Datasets Datasets/RazorAnalysis_SMCocktail_weighted_lumi-3.0_1-3btag_MultiJet.root

To make a .csv file of the yields in the cards directory, run the following command

    python python/WriteDataCard.py -b MultiJet -c config/run2_1-3btag.config -d cards Datasets/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_weighted_lumi-3.0_1-3btag_MultiJet.root Datasets/RazorInclusive_SMCocktail_weighted_lumi-3.0_1-3btag_MultiJet.root --print-yields 

To run 10000 "Bayesian" toys, marginalizing the shape parameters to
estimate the systematic uncertainty in each bin,

    python python/RunToys.py -b MultiJet -c config/run2_1-3btag.config -i FitResults/BinnedFitResults_MultiJet.root -d FitResults  -t 10000

More information on how to use the output from this toy generation
coming soon...
