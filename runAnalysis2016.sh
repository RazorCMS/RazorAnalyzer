#!/bin/sh
python python/ComputeScaleFactorsMacro.py 
python python/ComputeScaleFactorsNJetCorrectionByProcess.py 
python python/TestScaleFactors.py 
python python/TestGJetsScaleFactors.py 
root -l macros/BackgroundStudies/plotScaleFactorHistograms.C
python python/TTJetsDileptonCrossCheck.py
python python/VetoLeptonCrossCheck.py 
python python/DYJetsInvCrossCheck.py 
