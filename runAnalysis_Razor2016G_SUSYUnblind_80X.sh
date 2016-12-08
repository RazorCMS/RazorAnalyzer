#!/bin/sh
set -e
python python/ComputeScaleFactorsMacro.py --tag Razor2016G_SUSYUnblind_80X
python python/ComputeScaleFactorsNJetCorrectionByProcess.py --tag Razor2016G_SUSYUnblind_80X
python python/TestScaleFactors.py --tag Razor2016G_SUSYUnblind_80X
python python/TestGJetsScaleFactors.py --tag Razor2016G_SUSYUnblind_80X
python python/TTJetsDileptonCrossCheck.py --tag Razor2016G_SUSYUnblind_80X
python python/DYJetsInvCrossCheck.py --tag Razor2016G_SUSYUnblind_80X 
root -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C("Razor2016G_SUSYUnblind_80X")'
python python/DYJetsInvCrossCheck.py --tag Razor2016G_SUSYUnblind_80X --closure
python python/VetoLeptonCrossCheck.py --tag Razor2016G_SUSYUnblind_80X
python python/SignalRegionMacro.py --unblind --no-fit --tag Razor2016G_SUSYUnblind_80X
for nb in 0 1 2 3; do
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016G_SUSYUnblind_80X --box MultiJet --btags $nb
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016G_SUSYUnblind_80X --box LeptonMultiJet --btags $nb
done
for nb in 0 1 2; do
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016G_SUSYUnblind_80X --box DiJet --btags $nb
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016G_SUSYUnblind_80X --box LeptonJet --btags $nb
done
