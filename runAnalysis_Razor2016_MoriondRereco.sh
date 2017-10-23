#!/bin/sh
set -e
python python/ComputeScaleFactorsMacro.py 
python python/ComputeScaleFactorsNJetCorrectionByProcess.py 
python python/BTagClosureTestMacro.py
python python/TTJetsDileptonCrossCheck.py 
python python/GJetsBTagClosureTestMacro.py
python python/DYJetsInvCrossCheck.py 
root -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C("Razor2016_MoriondRereco")'
python python/DYJetsInvCrossCheck.py --closure
python python/VetoLeptonCrossCheck.py 
python python/SignalRegionMacro.py --unblind --fine-grained
python python/SignalRegionMacro.py --unblind --fine-grained --sideband
for nb in 0 1 2 3; do
    python python/SignalRegionPlotMacro.py --box MultiJet --btags $nb --fine-grained --unblind
    python python/SignalRegionPlotMacro.py --box LeptonMultiJet --btags $nb --fine-grained --unblind
    python python/SignalRegionPlotMacro.py --box MultiJet --btags $nb --fine-grained --unblind --sideband
    python python/SignalRegionPlotMacro.py --box LeptonMultiJet --btags $nb --fine-grained --unblind --sideband
done
for nb in 0 1 2; do
    python python/SignalRegionPlotMacro.py --box DiJet --btags $nb --fine-grained --unblind
    python python/SignalRegionPlotMacro.py --box LeptonJet --btags $nb --fine-grained  --unblind
    python python/SignalRegionPlotMacro.py --box DiJet --btags $nb --fine-grained --unblind --sideband
    python python/SignalRegionPlotMacro.py --box LeptonJet --btags $nb --fine-grained --unblind --sideband
done
