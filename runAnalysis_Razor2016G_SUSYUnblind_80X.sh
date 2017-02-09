#!/bin/sh
echo "===== ComputeScaleFactorsMacro.py ===="
python python/ComputeScaleFactorsMacro.py --tag Razor2016_MoriondRereco
echo "===== ComputeScaleFactorsNJetCorrectionByProcess.py ===="
python python/ComputeScaleFactorsNJetCorrectionByProcess.py --tag Razor2016_MoriondRereco
echo "===== TestScaleFactors.py ===="
python python/TestScaleFactors.py --tag Razor2016_MoriondRereco
echo "===== TestGJetsScaleFactors.py ===="
python python/TestGJetsScaleFactors.py --tag Razor2016_MoriondRereco
echo "===== TTJetsDileptonCrossCheck.py ===="
python python/TTJetsDileptonCrossCheck.py --tag Razor2016_MoriondRereco
echo "===== DYJetsInvCrossCheck.py ===="
python python/DYJetsInvCrossCheck.py --tag Razor2016_MoriondRereco 
echo "===== plotScaleFactorHistograms.C ===="
root -q -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C+g("Razor2016_MoriondRereco")'
echo "===== DYJetsInvCrossCheck.py --closure ===="
python python/DYJetsInvCrossCheck.py --tag Razor2016_MoriondRereco --closure
echo "===== VetoLeptonCrossCheck.py ===="
python python/VetoLeptonCrossCheck.py --tag Razor2016_MoriondRereco

echo "===== SIGNAL REGION ===="
python python/SignalRegionMacro.py --unblind --tag Razor2016_MoriondRereco
for nb in 0 1 2 3; do
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016_MoriondRereco --box MultiJet --btags $nb
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016_MoriondRereco --box LeptonMultiJet --btags $nb
done
for nb in 0 1 2; do
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016_MoriondRereco --box DiJet --btags $nb
    python python/SignalRegionPlotMacro.py --unblind --tag Razor2016_MoriondRereco --box LeptonJet --btags $nb
done
