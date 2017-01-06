#!/bin/bash

python python/DustinTuple2RooDataSet.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d Datasets_02Jan -b MultiJet --data -l 4400 ReReco2016_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_HTMHT_2016G_23Sep2016_SUSYUnblind.root

python python/DustinTuple2RooDataSet.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d Datasets_02Jan -b DiJet --data -l 4400 ReReco2016_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_HTMHT_2016G_23Sep2016_SUSYUnblind.root

python python/DustinTuple2RooDataSet.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d Datasets_02Jan -b LeptonMultiJet --data -l 4400 ReReco2016_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_SingleLepton_2016G_23Sep2016_SUSYUnblind.root

python python/DustinTuple2RooDataSet.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d Datasets_02Jan -b LeptonJet --data -l 4400 ReReco2016_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_SingleLepton_2016G_23Sep2016_SUSYUnblind.root

#MultiJet Full

python python/BinnedFit.py -b MultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_HTMHT_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-3btag_MultiJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/ -b MultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/BinnedFitResults_MultiJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/ -b MultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/BinnedFitResults_MultiJet.root  -t 1000

#python python/PlotFit.py -b MultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/BinnedFitResults_MultiJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/toys_Bayes_MultiJet.root -s fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/toys_Bayes_noStat_MultiJet.root --no-stat

python python/PlotFit.py -b MultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Full/BinnedFitResults_MultiJet.root --data -l 4400

#MultiJet Sideband

python python/BinnedFit.py -b MultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_HTMHT_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-3btag_MultiJet.root --data -l 4400 --fit-region LowMR,LowRsq 

#python python/BinnedFit.py -b MultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/ Datasets_02Jan/FullRazorInclusive_HTMHT_23Sep_23SepV1_Golden_lumi-4.400_0-3btag_MultiJet.root --data -l 4400 --fit-region LowMR,LowRsq 

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/ -b MultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/BinnedFitResults_MultiJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/ -b MultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/BinnedFitResults_MultiJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b MultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/BinnedFitResults_MultiJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/toys_Bayes_MultiJet.root -s fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/toys_Bayes_noStat_MultiJet.root --fit-region LowMR,LowRsq --no-stat

python python/PlotFit.py -b MultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/MultiJet/Sideband/BinnedFitResults_MultiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#DiJet Sideband

python python/BinnedFit.py -b DiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_HTMHT_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-2btag_DiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/ -b DiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/BinnedFitResults_DiJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/ -b DiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/BinnedFitResults_DiJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b DiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/BinnedFitResults_DiJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/toys_Bayes_DiJet.root -s fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/toys_Bayes_noStat_DiJet.root --fit-region LowMR,LowRsq --no-stat --print-errors

python python/PlotFit.py -b DiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Sideband/BinnedFitResults_DiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#DiJet Full

python python/BinnedFit.py -b DiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_HTMHT_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-2btag_DiJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/ -b DiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/BinnedFitResults_DiJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/ -b DiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/BinnedFitResults_DiJet.root  -t 1000

#python python/PlotFit.py -b DiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/BinnedFitResults_DiJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/toys_Bayes_DiJet.root -s fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/toys_Bayes_noStat_DiJet.root --no-stat --print-errors

python python/PlotFit.py -b DiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/DiJet/Full/BinnedFitResults_DiJet.root --data -l 4400

#LeptonMultiJet Sideband

python python/BinnedFit.py -b LeptonMultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_SingleLepton_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-3btag_LeptonMultiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/ -b LeptonMultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/BinnedFitResults_LeptonMultiJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/ -b LeptonMultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/BinnedFitResults_LeptonMultiJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b LeptonMultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/BinnedFitResults_LeptonMultiJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/toys_Bayes_LeptonMultiJet.root -s fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/toys_Bayes_noStat_LeptonMultiJet.root --fit-region LowMR,LowRsq --no-stat --print-errors

python python/PlotFit.py -b LeptonMultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Sideband/BinnedFitResults_LeptonMultiJet.root --data -l 4400 --fit-region LowMR,LowRsq

#LeptonMultiJet Full

python python/BinnedFit.py -b LeptonMultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_SingleLepton_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-3btag_LeptonMultiJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/ -b LeptonMultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/BinnedFitResults_LeptonMultiJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/ -b LeptonMultiJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/BinnedFitResults_LeptonMultiJet.root  -t 1000

#python python/PlotFit.py -b LeptonMultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/BinnedFitResults_LeptonMultiJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/toys_Bayes_LeptonMultiJet.root -s fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/toys_Bayes_noStat_LeptonMultiJet.root --no-stat --print-errors

python python/PlotFit.py -b LeptonMultiJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonMultiJet/Full/BinnedFitResults_LeptonMultiJet.root --data -l 4400

#LeptonJet Sideband

python python/BinnedFit.py -b LeptonJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_SingleLepton_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-2btag_LeptonJet.root --data -l 4400 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/ -b LeptonJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/BinnedFitResults_LeptonJet.root --no-stat -t 1000 --fit-region LowMR,LowRsq

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/ -b LeptonJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/BinnedFitResults_LeptonJet.root  -t 1000 --fit-region LowMR,LowRsq

#python python/PlotFit.py -b LeptonJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/BinnedFitResults_LeptonJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/toys_Bayes_LeptonJet.root -s fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/toys_Bayes_noStat_LeptonJet.root --fit-region LowMR,LowRsq --no-stat --print-errors

python python/PlotFit.py -b LeptonJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Sideband/BinnedFitResults_LeptonJet.root --data -l 4400 --fit-region LowMR,LowRsq

#LeptonJet Full

python python/BinnedFit.py -b LeptonJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/ Datasets_02Jan/FullRazorInclusive_Razor2016_MoriondRereco_SingleLepton_2016G_23Sep2016_SUSYUnblind_lumi-4.400_0-2btag_LeptonJet.root --data -l 4400 

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/ -b LeptonJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/BinnedFitResults_LeptonJet.root --no-stat -t 1000

#python python/RunToys.py -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/ -b LeptonJet -l 4400 -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/BinnedFitResults_LeptonJet.root  -t 1000

#python python/PlotFit.py -b LeptonJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/BinnedFitResults_LeptonJet.root --data -l 4400 -t fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/toys_Bayes_LeptonJet.root -s fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/toys_Bayes_noStat_LeptonJet.root --no-stat --print-errors

python python/PlotFit.py -b LeptonJet -c config/run2_2017_01_03_Run2016G_SUSYUnblind_Sep23ReReco.config -d fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/ -i fits_2017_01_03/ReReco2016_02Jan/LeptonJet/Full/BinnedFitResults_LeptonJet.root --data -l 4400
