#!/bin/sh
set -e
python python/ComputeScaleFactorsMacro.py --tag Razor2016G_SUSYUnblind_80X
python python/ComputeScaleFactorsNJetCorrectionByProcess.py --tag Razor2016G_SUSYUnblind_80X
python python/TestScaleFactors.py --tag Razor2016G_SUSYUnblind_80X
python python/TestGJetsScaleFactors.py --tag Razor2016G_SUSYUnblind_80X
python python/TTJetsDileptonCrossCheck.py --tag Razor2016G_SUSYUnblind_80X
python python/DYJetsInvCrossCheck.py --tag Razor2016G_SUSYUnblind_80X
root -l 'macros/BackgroundStudies/plotScaleFactorHistograms.C("Razor2016G_SUSYUnblind_80X")'
python python/VetoLeptonCrossCheck.py --tag Razor2016G_SUSYUnblind_80X
