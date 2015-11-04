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

    config/run2_sideband.config

Setup combine from lxplus:

	mkdir ~/work/RAZORRUN2/
	cd ~/work/RAZORRUN2/
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

    mkdir Datasets
    mkdir FitResults; mkdir FitProjections; mkdir cards

Now convert the SM MC ntuples into a SM Cocktail RooDataSet (ignoring
the QCD contirbution).

     python python/DustinTuple2RooDataSet.py -c config/run2_sideband.config -b MultiJet -d Datasets/ -w -l 4000 \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_DYJetsToLL_M-50_HTBinned_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_DYJetsToLL_M-5to50_HTBinned_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_SingleTop_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_TTV_1pb_weighted.rootRazorInclusive_VV_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_WJetsToLNu_HTBinned_1pb_weighted.root \
	 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p19_ForFullStatus20151030/MC/RazorInclusive_ZJetsToNuNu_HTBinned_1pb_weighted.root 

Similarly you can run over the data,

     python python/DustinTuple2RooDataSet.py -b MultiJet -c config/run2_sideband.config -d Datasets/ --data -l 1264 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/RazorInclusive/V1p20_ForFullStatus20151030/Data/RazorInclusive_HTMHT_Run2015D_Oct05ReMiniAOD_PRv4_GoodLumiGolden.root

To perform the fit on the SM Cocktail,

     python python/BinnedFit.py -b MultiJet -c config/run2_sideband.config -d FitResults Datasets/RazorInclusive_SMCocktail_weighted_lumi-4.000_0-3btag_MultiJet.root

To produce the signal templates,

	 python python/SMSTemplates.py -c config/run2_sideband.config -b MultiJet -d Datasets/ -l 4000 root://eoscms.cern.ch//eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/V1p22_ForPreappFreezing20151106/jobs/combined/SMS-T1bbbb_1500_100.root 

To run Bayesian toys, marginalizing the shape parameters to
estimate the systematic uncertainty in each bin, and produce a large number of other 1D fit projections, and the 2D
projections in MR and Rsq,

     python python/RunToys.py -b MultiJet -c config/run2_sideband.config -i FitResults/BinnedFitResults_MultiJet.root -d FitResults  -t 10000
     python python/PlotFit.py -b MultiJet -c config/run2_sideband.config -d FitResults -i FitResults/BinnedFitResults_MultiJet.root -t FitResults/toys_Bayes_MultiJet.root
	
Next, to produce the datacards and run combine, execute the following:

     python python/DustinTuple2RooDataSet.py -b MultiJet -c config/run2_sideband.config -w -l 4000 -d Datasets Signals/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_1pb_weighted.root
     python python/RunCombine.py -b MultiJet -c config/run2_sideband.config -d cards --lumi-array 4 -m T1bbbb --mGluino 1500 --mLSP 100

To submit jobs to run the limits on the data:

	 mkdir ~/work/RAZORRUN2/Limits/
     python python/RunCombineJobs.py -i FitResults/BinnedFitResults.root -c config/run2_sideband.config --data -l 1.264 -m T1bbbb -b MultiJet -d cards/ -q 8nh

After the jobs finish, the output files will be in Limits directory we
just created

	 cp ~/work/RAZORRUN2/Limits/cards/higgs*.root cards/
     python python/GetCombine.py -m T1bbbb -b MultiJet -c config/run2_sideband.config -d cards/ -l 1.264
     python python/Get2DContour.py -m T1bbbb -b MultiJet  -d cards/ 

To make the final (expected) limit plot, we need to check out a different repository. 

     git clone git@github.com:RazorCMS/PlotsSMS
     cd PlotsSMS

Note you have to change the smoothed cross section limit file location
in config/SUS15004/T1bbbb\_Exp\_SUS15004.cfg. Then you can run it,

	 python python/makeSMSPlots.py config/SUS15004/T1bbbb_Exp_SUS15004.cfg T1bbbbAsymptotic

To make an "unweighted" dataset (not needed for the preceding commands),
execute

     python python/RooDataSet2UnweightedDataSet.py -b MultiJet -c config/run2_sideband.config -d Datasets Datasets/RazorAnalysis_SMCocktail_weighted_lumi-4.000_0-3btag_MultiJet.root

To make a .csv file of the yields in the cards directory, run the following command

     python python/WriteDataCard.py -b MultiJet -c config/run2_sideband.config -d cards Datasets/RazorInclusive_SMS-T1bbbb_2J_mGl-1500_mLSP-100_weighted_lumi-4.000_0-3btag_MultiJet.root Datasets/RazorInclusive_SMCocktail_weighted_lumi-3.0_1-3btag_MultiJet.root --print-yields 

