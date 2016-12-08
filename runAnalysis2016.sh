#!/bin/sh
set -e 
python python/ComputeScaleFactorsMacro.py 
python python/ComputeScaleFactorsNJetCorrectionByProcess.py 
python python/TestScaleFactors.py 
python python/TestGJetsScaleFactors.py 
python python/TTJetsDileptonCrossCheck.py 
python python/DYJetsInvCrossCheck.py 
root -l macros/BackgroundStudies/plotScaleFactorHistograms.C
python python/DYJetsInvCrossCheck.py --closure
python python/VetoLeptonCrossCheck.py 
python python/SignalRegionMacro.py --unblind --no-fit 
for nb in 0 1 2 3; do
    python python/SignalRegionPlotMacro.py --unblind --box MultiJet --btags $nb
    python python/SignalRegionPlotMacro.py --unblind --box LeptonMultiJet --btags $nb
done
for nb in 0 1 2; do
    python python/SignalRegionPlotMacro.py --unblind --box DiJet --btags $nb
    python python/SignalRegionPlotMacro.py --unblind --box LeptonJet --btags $nb
done
