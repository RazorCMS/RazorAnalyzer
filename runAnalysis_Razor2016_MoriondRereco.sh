#!/bin/sh
set -e
python python/ComputeScaleFactorsMacro.py --tag Razor2016_MoriondRereco
python python/ComputeScaleFactorsNJetCorrectionByProcess.py --tag Razor2016_MoriondRereco
python python/TestScaleFactors.py --tag Razor2016_MoriondRereco
python python/TestGJetsScaleFactors.py --tag Razor2016_MoriondRereco
python python/TTJetsDileptonCrossCheck.py --tag Razor2016_MoriondRereco
python python/DYJetsInvCrossCheck.py --tag Razor2016_MoriondRereco 
root -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C("Razor2016_MoriondRereco")'
python python/DYJetsInvCrossCheck.py --tag Razor2016_MoriondRereco --closure
python python/VetoLeptonCrossCheck.py --tag Razor2016_MoriondRereco
python python/SignalRegionMacro.py --unblind --tag Razor2016_MoriondRereco --no-fit
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
