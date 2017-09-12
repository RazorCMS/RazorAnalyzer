#!/bin/sh
set -e
#python python/ComputeScaleFactorsMacroDM.py --tag Razor2016_MoriondRereco
#python python/ComputeScaleFactorsNJetCorrectionByProcessDM.py --tag Razor2016_MoriondRereco
#python python/TestScaleFactorsDM.py --tag Razor2016_MoriondRereco
#echo "===== TestGJetsScaleFactorsDM.py ===="
#python python/TestGJetsScaleFactorsDM.py --tag Razor2016_MoriondRereco
#python python/TTJetsDileptonCrossCheckDM.py --tag Razor2016_MoriondRereco
#python python/DYJetsInvCrossCheckDM.py --tag Razor2016_MoriondRereco 
#root -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C("Razor2016_MoriondRereco")'
#python python/DYJetsInvCrossCheckDM.py --tag Razor2016_MoriondRereco --closure
#python python/VetoLeptonCrossCheckDM.py --tag Razor2016_MoriondRereco
#python python/TestGJetsScaleFactorsDM.py --tag Razor2016_MoriondRereco
python python/SignalRegionMacroDM.py --tag Razor2016_MoriondRereco --no-fit --verbose
#for box in MultiJet DiJet LeptonJet LeptonMultiJet; do
#    python python/AddRazorFitToMADDResults.py --tag Razor2016_MoriondRereco --box $box
#done
python python/AddRazorFitToMADDResultsDM.py --tag Razor2016_MoriondRereco --box MultiJet
python python/SignalRegionPlotMacroDM.py --tag Razor2016_MoriondRereco --box MultiJet --btags 0
#for nb in 0 1 2 3; do
#    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box MultiJet --btags $nb
#    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box LeptonMultiJet --btags $nb
#done
#for nb in 0 1 2; do
#    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box DiJet --btags $nb
#    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box LeptonJet --btags $nb
#done
