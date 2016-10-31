#!/bin/sh
python python/ComputeScaleFactorsMacro.py 
python python/ComputeScaleFactorsNJetCorrection.py 
python python/TestScaleFactors.py 
python python/TestGJetsScaleFactors.py 
root -l macros/BackgroundStudies/plotScaleFactorHistograms.C
root -l macros/BackgroundStudies/TTBarDileptonCrossCheck.C++(1)
root -l macros/BackgroundStudies/TTBarDileptonCrossCheck.C++(2)
python python/VetoLeptonCrossCheck.py 
python python/DYJetsInvCrossCheck.py 
