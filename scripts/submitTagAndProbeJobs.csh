#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================



##########################
# Tag And Probe
##########################
foreach sample( \
SingleElectron_Run2015B \
DoubleEG_Run2015B \
SingleMuon_Run2015B \
DoubleMuon_Run2015B \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/data/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10551 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27LooseEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10552 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27TightEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10553 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle32TightEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10560 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 10561 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedExtEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20551 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu20EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20552 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20553 $filesPerJob $jobnumber TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerMu50EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist true 20560 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27ORMu50EffDenominatorTight/jobs/


  end
end

foreach sample( \
DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8\
)
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/Run2/razorNtuplerV1p16/MC/50ns/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10551 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27LooseEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27LooseEffDenominatorTight/jobs/ 
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10552 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle27TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle27TightEffDenominatorTight/jobs/ 
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10553 $filesPerJob $jobnumber TagAndProbe_EleTriggerEle32TightEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEle32TightEffDenominatorTight/jobs/ 
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10560 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedEffDenominatorTight/jobs/
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 10561 $filesPerJob $jobnumber TagAndProbe_EleTriggerEleCombinedExtEffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/EleTriggerEleCombinedExtEffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20551 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu20EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20552 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu27EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20553 $filesPerJob $jobnumber TagAndProbe_MuTriggerMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerMu50EffDenominatorTight/jobs/

	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber}.out -J TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RazorTagAndProbe $inputfilelist false 20560 $filesPerJob $jobnumber TagAndProbe_MuTriggerIsoMu20ORMu50EffDenominatorTight_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TagAndProbe/MuTriggerIsoMu27ORMu50EffDenominatorTight/jobs/


  end
end


