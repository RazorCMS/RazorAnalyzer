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
python python/SignalRegionMacro.py --unblind --fine-grained
for nb in 0 1 2 3; do
    python python/SignalRegionPlotMacro.py --box MultiJet --btags $nb --fine-grained 
    python python/SignalRegionPlotMacro.py --box LeptonMultiJet --btags $nb --fine-grained 
done
for nb in 0 1 2; do
    python python/SignalRegionPlotMacro.py --box DiJet --btags $nb --fine-grained 
    python python/SignalRegionPlotMacro.py --box LeptonJet --btags $nb --fine-grained 
done
