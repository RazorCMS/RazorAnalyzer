#!/bin/sh
set -e
python python/ComputeScaleFactorsMacro.py 
python python/ComputeScaleFactorsNJetCorrectionByProcess.py 
python python/TestScaleFactors.py 
python python/TTJetsDileptonCrossCheck.py 
python python/DYJetsInvCrossCheck.py 
root -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C("Razor2016_MoriondRereco")'
python python/DYJetsInvCrossCheck.py --closure
python python/VetoLeptonCrossCheck.py 
python python/TestGJetsScaleFactors.py 
python python/SignalRegionMacro.py --unblind --no-fit
for box in MultiJet DiJet LeptonJet LeptonMultiJet; do
    python python/AddRazorFitToMADDResults.py --tag Razor2016_MoriondRereco --box $box
done
for nb in 0 1 2 3; do
    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box MultiJet --btags $nb
    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box LeptonMultiJet --btags $nb
done
for nb in 0 1 2; do
    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box DiJet --btags $nb
    python python/SignalRegionPlotMacro.py --tag Razor2016_MoriondRereco --box LeptonJet --btags $nb
done
