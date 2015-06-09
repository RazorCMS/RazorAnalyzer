#!/bin/tcsh
#===================================================================================================
# Submit a set of jobs to run over a given dataset, splitting the jobs according to the filesets.
#
# Version 1.0                                                                      November 14, 2008
#===================================================================================================


##########################
# Baseline Razor Analysis
##########################
foreach sample( \
SMS-T1tttt_2J_mGl-1500_mLSP-100_20bx25 \
SMS-T1bbbb_2J_mGl-1500_mLSP-100_20bx25 \
SMS-T1tttt_2J_mGl-1200_mLSP-800_20bx25 \
SMS-T1bbbb_2J_mGl-1000_mLSP-900_20bx25 \
SMS-T1qqqq_2J_mGl-1400_mLSP-100_20bx25 \
SMS-T1qqqq_2J_mGl-1000_mLSP-800_20bx25 \
SMS-T2bb_2J_mStop-600_mLSP-580_20bx25 \
SMS-T2bb_2J_mStop-900_mLSP-100_20bx25 \
SMS-T2qq_2J_mStop-1200_mLSP-100_20bx25 \
SMS-T2qq_2J_mStop-600_mLSP-550_20bx25 \
SMS-T2tt_2J_mStop-425_mLSP-325_20bx25 \
SMS-T2tt_2J_mStop-500_mLSP-325_20bx25 \
SMS-T2tt_2J_mStop-650_mLSP-325_20bx25 \
SMS-T2tt_2J_mStop-850_mLSP-100_20bx25 \
TTJets_20bx25 \
QCD_HT100To250_20bx25 \
QCD_HT250To500_20bx25 \
QCD_HT500To1000_20bx25 \
QCD_HT1000ToInf_20bx25 \
WJetsToLNu_HT100To200_20bx25 \
WJetsToLNu_HT200To400_20bx25 \
WJetsToLNu_HT400To600_20bx25 \
WJetsToLNu_HT600ToInf_20bx25 \
DYJetsToLL_HT100To200_20bx25 \
DYJetsToLL_HT200To400_20bx25 \
DYJetsToLL_HT400To600_20bx25 \
DYJetsToLL_HT600ToInf_20bx25 \
ZJetsToNuNu_HT100To200_20bx25 \
ZJetsToNuNu_HT200To400_20bx25 \
ZJetsToNuNu_HT400To600_20bx25 \
ZJetsToNuNu_HT600ToInf_20bx25 \
T_tW_20bx25 \
Tbar_tW_20bx25 \
TToLeptons_s_20bx25 \
TBarToLeptons_s_20bx25 \
TToLeptons_t_20bx25 \
TBarToLeptons_t_20bx25 \
WZJetsTo3LNu_20bx25 \
TTWJets_20bx25 \
TTZJets_20bx25 \
SMS-T1qqqq_2J_mGl-1000_mLSP-800_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14 \
SMS-T1bbbb_2J_mGl-1000_mLSP-900_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
SMS-T2tt_2J_mStop-425_mLSP-325_Tune4C_13TeV-madgraph-tauola_AVE30BX50_ST_PHYS14
DYJetsToLL_M-50_13TeV-madgraph-pythia8-tauola_AVE30BX50_ST_PHYS14
DYJetsToLL_M-50_13TeV-madgraph-pythia8_4bx50_PHYS14
TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_30bx50_PHYS14
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p5-25ns_v1_v7/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorInclusive_${jobnumber}.out -J RazorInclusive_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh razor $inputfilelist false -1 $filesPerJob $jobnumber RazorInclusive_${sample}_20bx25.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RazorInclusive/jobs/
    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HbbRazor_${jobnumber}.out -J HbbRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hbbrazor $inputfilelist false -1 $filesPerJob $jobnumber HbbRazor_${sample}_20bx25.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HbbRazor/jobs/
    sleep 0.1
  end

end

##########################
# Run One RazorInclusive
##########################
foreach sample( \
Data_DoubleElectron_Run2012A \
Data_DoubleElectron_Run2012B \
Data_DoubleElectron_Run2012C \
Data_DoubleElectron_Run2012D \
Data_DoubleMuParked_Run2012A \
Data_DoubleMuParked_Run2012B \
Data_DoubleMuParked_Run2012C \
Data_DoubleMuParked_Run2012D \
Data_MuEG_Run2012A \
Data_MuEG_Run2012B \
Data_MuEG_Run2012C \
Data_MuEG_Run2012D \
Data_SingleElectron_Run2012A \
Data_SingleElectron_Run2012B \
Data_SingleElectron_Run2012C \
Data_SingleElectron_Run2012D \
Data_SingleMu_Run2012A \
Data_SingleMu_Run2012B \
Data_SingleMu_Run2012C \
Data_SingleMu_Run2012D \
Data_HT_Run2012A \
Data_HTMHTParked_Run2012B \
Data_HTMHTParked_Run2012C \
Data_HTMHTParked_Run2012D \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p8-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run2/RunOneRazorInclusive/jobs/RazorInclusive_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/public/condor/res/Run2SUSY/RazorAnalysis/RazorInclusive_${sample}_${jobnumber}.out -J RazorInclusive_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazor $inputfilelist true 1 $filesPerJob $jobnumber RazorInclusive_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/run2/RunOneRazorInclusive/jobs/
	
        while ( `bjobs | grep RUN | wc | awk '{print $1}' | awk '{if($1 > 50){print "yes"}else{print "no"}}'` == "yes")
	   sleep 10
        end
    endif
  end
end

foreach sample( \
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball \
WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball \
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola \
T_s-channel_TuneZ2star_8TeV-powheg-tauola \
T_t-channel_TuneZ2star_8TeV-powheg-tauola \
T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo3LNu_8TeV_TuneZ2Star_madgraph_tauola \
ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola \
TTWJets_8TeV-madgraph \
TTWWJets_8TeV-madgraph \
TTZJets_8TeV-madgraph_v2 \
TTTT_TuneZ2star_8TeV-madgraph-tauola \
WJetsToLNu_HT-150To200_8TeV-madgraph \
WJetsToLNu_HT-200To250_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph_v2 \
WJetsToLNu_HT-300To400_8TeV-madgraph \
WJetsToLNu_HT-300To400_8TeV-madgraph_v2 \
WJetsToLNu_HT-400ToInf_8TeV-madgraph \
WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2 \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
TTJets_HadronicMGDecays_8TeV-madgraph \
TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola \
QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia \
QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6 \
QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6 \
QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6 \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_ext \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p8-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/run2/RunOneRazorInclusive/jobs/RazorInclusive_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorInclusive_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazor $inputfilelist false 1 $filesPerJob $jobnumber RazorInclusive_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/run2/RunOneRazorInclusive/jobs/
	while ( `bjobs | grep RUN | wc | awk '{print $1}' | awk '{if($1 > 50){print "yes"}else{print "no"}}'` == "yes")
	   sleep 10
        end
    endif
  end
end



##########################
# Veto Lepton Study
##########################
foreach sample( \
T1bbbb_1500
T1tttt_1500
TTJets \
) 
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/lists/razorNtuplerV1p4-25ns_v1_v1/${sample}_20bx25.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/VetoLeptonStudy_${jobnumber}.out -J RazorAnalysis_VetoLeptonStudy_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorVetoLeptonStudy $inputfilelist 1 $filesPerJob $jobnumber VetoLeptonStudy_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/VetoLeptonStudy/jobs/
    sleep 0.1
  end

end


##########################
# Object Ntuplers
##########################
foreach sample( \
TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Spring15_25ns \
)
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p10-Run2/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/ElectronNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/ElectronNtupler_${jobnumber}.out -J RazorAnalysis_ElectronNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh electronNtupler $inputfilelist false 1 $filesPerJob $jobnumber ElectronNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/  
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/ElectronNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/ElectronNtupler_${jobnumber}.out -J RazorAnalysis_ElectronNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh electronNtupler $inputfilelist false 0 $filesPerJob $jobnumber ElectronNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/ElectronNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/ElectronNtupler_${jobnumber}.out -J RazorAnalysis_ElectronNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh electronNtupler $inputfilelist false 11 $filesPerJob $jobnumber ElectronNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/ElectronNtuple/jobs/
      sleep 0.1
    endif 
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/MuonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/MuonNtupler_${jobnumber}.out -J RazorAnalysis_MuonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh muonNtupler $inputfilelist false 1 $filesPerJob $jobnumber MuonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/MuonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/MuonNtupler_${jobnumber}.out -J RazorAnalysis_MuonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh muonNtupler $inputfilelist false 0 $filesPerJob $jobnumber MuonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/MuonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/MuonNtupler_${jobnumber}.out -J RazorAnalysis_MuonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh muonNtupler $inputfilelist false 11 $filesPerJob $jobnumber MuonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/MuonNtuple/jobs/
      sleep 0.1
    endif 

    bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/JetNtupler_${jobnumber}.out -J RazorAnalysis_JetNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh jetNtupler $inputfilelist false -1 $filesPerJob $jobnumber JetNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/JetNtuple/jobs/
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TauNtupler_${jobnumber}.out -J RazorAnalysis_TauNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh tauNtupler $inputfilelist false 1 $filesPerJob $jobnumber TauNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TauNtuple/jobs/
    sleep 0.1
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TauNtupler_${jobnumber}.out -J RazorAnalysis_TauNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh tauNtupler $inputfilelist false 0 $filesPerJob $jobnumber TauNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TauNtuple/jobs/
    sleep 0.1
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/TauNtupler_${jobnumber}.out -J RazorAnalysis_TauNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh tauNtupler $inputfilelist false 11 $filesPerJob $jobnumber TauNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/TauNtuple/jobs/
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/PhotonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/PhotonNtupler_${jobnumber}.out -J RazorAnalysis_PhotonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh photonNtupler $inputfilelist false 1 $filesPerJob $jobnumber PhotonNtuple_Prompt_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/PhotonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/PhotonNtupler_${jobnumber}.out -J RazorAnalysis_PhotonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh photonNtupler $inputfilelist false 0 $filesPerJob $jobnumber PhotonNtuple_Fake_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/
      sleep 0.1
    endif 

    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/PhotonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root ) then
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/PhotonNtupler_${jobnumber}.out -J RazorAnalysis_PhotonNtupler_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh photonNtupler $inputfilelist false 11 $filesPerJob $jobnumber PhotonNtuple_PromptGenLevel_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/PhotonNtuple/jobs/
      sleep 0.1
    endif 

  end
end



##########################
# Control Region Study
##########################
foreach sample( \
SMS-T1tttt_2J_mGl-1500_mLSP-100\
SMS-T1bbbb_2J_mGl-1500_mLSP-100\
SMS-T1tttt_2J_mGl-1200_mLSP-800\
SMS-T1bbbb_2J_mGl-1000_mLSP-900\
SMS-T1qqqq_2J_mGl-1400_mLSP-100\
SMS-T1qqqq_2J_mGl-1000_mLSP-800\
TTJets \
WJetsToLNu_HT100To200 \
WJetsToLNu_HT200To400 \
WJetsToLNu_HT400To600 \
WJetsToLNu_HT600ToInf \
DYJetsToLL_HT100To200 \
DYJetsToLL_HT200To400 \
DYJetsToLL_HT400To600 \
DYJetsToLL_HT600ToInf \
T_tW \
TBar_tW \
TToLeptons_s \
TBarToLeptons_s \
TToLeptons_t \
TBarToLeptons_t \
WZJetsTo3LNu \
TTWJets \
TTZJets \
QCDHT100To250 \
QCDHT250To500 \
QCDHT500To1000 \
QCDHT1000ToInf \
) 
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/lists/razorNtuplerV1p4-25ns_v1_v1/${sample}_20bx25.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh RazorControlRegions $inputfilelist 4 $filesPerJob $jobnumber RazorControlRegions_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorControlRegions/jobs/
    sleep 0.1
  end
  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/VetoLeptonEfficiencyDileptonControlRegion_${jobnumber}.out -J RazorAnalysis_VetoLeptonEfficiencyDileptonControlRegion_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh VetoLeptonEfficiencyControlRegion $inputfilelist 0 $filesPerJob $jobnumber VetoLeptonEfficiencyDileptonControlRegion_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/VetoLeptonEfficiencyDileptonControlRegion/jobs/
    sleep 0.1
  end
  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/VetoLeptonEfficiencySingleLeptonControlRegion_${jobnumber}.out -J RazorAnalysis_VetoLeptonEfficiencySingleLeptonControlRegion_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_2_0/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh VetoLeptonEfficiencyControlRegion $inputfilelist 1 $filesPerJob $jobnumber VetoLeptonEfficiencySingleLeptonControlRegion_${sample}_25ns.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/VetoLeptonEfficiencySingleLeptonControlRegion/jobs/
    sleep 0.1
  end

end



##########################
# Run One Control Region Study
##########################
foreach sample( \
SingleMu \
SingleElectron \
DoubleMuParked \
DoubleElectron \
MuEG \
HT \
HTMHTParked \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazorControlRegions $inputfilelist true 1 $filesPerJob $jobnumber RunOneRazorControlRegions_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazorControlRegions $inputfilelist true 2 $filesPerJob $jobnumber RunOneRazorControlRegions_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end


 foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_RazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh RunOneRazorControlRegions $inputfilelist true 10 $filesPerJob $jobnumber RunOneRazorControlRegions_RazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_SingleLeptonRazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazorControlRegions $inputfilelist true 12 $filesPerJob $jobnumber RunOneRazorControlRegions_SingleLeptonRazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end

end



foreach sample( \
DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball \
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola \
TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola \
TTJets_HadronicMGDecays_8TeV-madgraph \
WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball \
WJetsToLNu_HT-150To200_8TeV-madgraph \
WJetsToLNu_HT-200To250_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph_v2 \
WJetsToLNu_HT-300To400_8TeV-madgraph_v2 \
WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2 \
T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
T_s-channel_TuneZ2star_8TeV-powheg-tauola \
T_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola\
WZJetsTo3LNu_8TeV_TuneZ2Star_madgraph_tauola \
ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
TTWJets_8TeV-madgraph \
TTWWJets_8TeV-madgraph \
TTZJets_8TeV-madgraph_v2 \
TTTT_TuneZ2star_8TeV-madgraph-tauola \
) 
  setenv LSB_JOB_REPORT_MAIL N
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazorControlRegions $inputfilelist false 1 $filesPerJob $jobnumber RunOneRazorControlRegions_DileptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazorControlRegions $inputfilelist false 2 $filesPerJob $jobnumber RunOneRazorControlRegions_SingleLeptonSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end


 foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_RazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh RunOneRazorControlRegions $inputfilelist false 10 $filesPerJob $jobnumber RunOneRazorControlRegions_RazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/RunOneRazorControlRegions_SingleLeptonRazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorControlRegions_${sample}_${jobnumber}.out -J RazorAnalysis_RazorControlRegions_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh RunOneRazorControlRegions $inputfilelist false 12 $filesPerJob $jobnumber RunOneRazorControlRegions_SingleLeptonRazorSkim_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/RunOneRazorControlRegions/jobs/
	sleep 0.1
    endif
  end

end







##########################
# Razor Z Analysis
##########################
foreach sample( \
#DYJetsToLL_MG \
#DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYToEE_powheg
DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
Data_DoubleMuParked_Run2012A \
Data_DoubleMuParked_Run2012B \
Data_DoubleMuParked_Run2012C \
Data_DoubleMuParked_Run2012D \
Data_DoubleElectron_Run2012A \
Data_DoubleElectron_Run2012B \
Data_DoubleElectron_Run2012C \
Data_DoubleElectron_Run2012D \
DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6 \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball \
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola \
T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
T_s-channel_TuneZ2star_8TeV-powheg-tauola \
T_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo3LNu_8TeV_TuneZ2Star_madgraph_tauola \
ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
TTWJets_8TeV-madgraph \
TTWWJets_8TeV-madgraph \
TTZJets_8TeV-madgraph_v2 \
TTTT_TuneZ2star_8TeV-madgraph-tauola \
) 
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/RazorZAnalysis_${sample}_${jobnumber}.out -J RazorAnalysis_RazorZAnalysis_${sample}_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorZAnalysis $inputfilelist 1 $filesPerJob $jobnumber RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/
	sleep 0.1
    endif
  end

end




  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh RazorZAnalysis $inputfilelist 1 $filesPerJob $jobnumber RazorZAnalysis_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/
	sleep 0.1
    endif
  end

 foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/RazorZAnalysis_DielectronSkim_${sample}.Job${jobnumber}Of${maxjob}.root ) then
	echo "job " $jobnumber " out of " $maxjob
	/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh RazorZAnalysis $inputfilelist 2 $filesPerJob $jobnumber RazorZAnalysis_DielectronSkim_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/user/s/sixie/work/public/Run2SUSY/RazorZAnalysis/jobs/
	sleep 0.1
    endif
  end





##########################
# HZZ Razor Analyzer
##########################
foreach sample( \
ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball \
TTJets_FullLeptMGDecays_8TeV-madgraph-tauola \
T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
T_s-channel_TuneZ2star_8TeV-powheg-tauola \
T_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola \
Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola \
WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola\
WZJetsTo3LNu_8TeV_TuneZ2Star_madgraph_tauola \
ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola \
ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola \
TTWJets_8TeV-madgraph \
TTWWJets_8TeV-madgraph \
TTZJets_8TeV-madgraph_v2 \
TTTT_TuneZ2star_8TeV-madgraph-tauola \
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p8-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root ) then
      echo "job " $jobnumber " out of " $maxjob
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HZZRazor_${jobnumber}.out -J HZZRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hzzRunOneRazor $inputfilelist false -1 $filesPerJob $jobnumber HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/    
    endif 
    sleep 0.1
  end

end

foreach sample( \
GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6 \
TTbarH_HToZZTo4L_M-125_8TeV-pythia6 \
VBF_HToZZTo4L_M-125_8TeV-powheg-pythia6 \
WH_HToZZTo4L_M-125_8TeV-pythia6 \
ZH_HToZZTo4L_M-125_8TeV-pythia6 \
GluGluToZZTo4L_8TeV-gg2zz-pythia6 \
GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6 \
ZZTo2e2mu_8TeV-powheg-pythia6 \
ZZTo2e2mu_8TeV_ext-powheg-pythia6 \
ZZTo2e2tau_8TeV-powheg-pythia6 \
ZZTo2e2tau_8TeV_ext-powheg-pythia6 \
ZZTo2mu2tau_8TeV-powheg-pythia6 \
ZZTo2mu2tau_8TeV_ext-powheg-pythia6 \
ZZTo4e_8TeV-powheg-pythia6 \
ZZTo4e_8TeV_ext-powheg-pythia6 \
ZZTo4mu_8TeV-powheg-pythia6 \
ZZTo4mu_8TeV_ext-powheg-pythia6 \
ZZTo4tau_8TeV-powheg-pythia6 \
ZZTo4tau_8TeV_ext-powheg-pythia6 \
GluGluToHToZZTo4L_M-125_8TeV-powheg15-pythia6 \
ZZTo4e_8TeV_mll8_mZZ95-160-powheg15-pythia6 \
ZZTo4mu_8TeV_mll8_mZZ95-160-powheg15-pythia6 \
ZZTo4tau_8TeV_mll8_mZZ95-160-powheg15-pythia6 \
GluGluToHToZZTo4L_M-125_8TeV-minloHJJ-pythia6-tauola \
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p9-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root ) then
      echo "job " $jobnumber " out of " $maxjob
      bsub -q 1nd -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HZZRazor_${jobnumber}.out -J HZZRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hzzRunOneRazor $inputfilelist false -1 $filesPerJob $jobnumber HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/    
    endif
    sleep 0.1
  end

end


foreach sample( \
Data_DoubleMuParked_Run2012A \
Data_DoubleMuParked_Run2012B \
Data_DoubleMuParked_Run2012C \
Data_DoubleMuParked_Run2012D \
Data_DoubleElectron_Run2012A \
Data_DoubleElectron_Run2012B \
Data_DoubleElectron_Run2012C \
Data_DoubleElectron_Run2012D \
Data_MuEG_Run2012A \
Data_MuEG_Run2012B \
Data_MuEG_Run2012C \
Data_MuEG_Run2012D \
) 
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/lists/razorNtuplerV1p8-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`

  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root ) then
      echo "job " $jobnumber " out of " $maxjob
      bsub -q 8nm -o /afs/cern.ch/user/s/sixie/work/private/condor/res/Run2SUSY/RazorAnalysis/HZZRazor_${jobnumber}.out -J HZZRazor_Razor_${jobnumber} /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh hzzRunOneRazor $inputfilelist true -1 $filesPerJob $jobnumber HZZRazor_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/HZZRazor/jobs/    
      sleep 0.1
    endif
  end

end






##########################
# MC for ZInv
##########################
foreach sample( \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
WJetsToLNu_HT-150To200_8TeV-madgraph \
WJetsToLNu_HT-200To250_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph_v2 \
WJetsToLNu_HT-300To400_8TeV-madgraph \
WJetsToLNu_HT-300To400_8TeV-madgraph_v2 \
WJetsToLNu_HT-400ToInf_8TeV-madgraph \
WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2 \
) 
  set inputfilelist="/afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/RazorControlRegions_${jobnumber}.out -J RazorAnalysis_RazorPhoton_${jobnumber} /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorPhotonStudy $inputfilelist 0 $filesPerJob $jobnumber RazorControlRegions_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/jobs/
    sleep 0.1
  end

end

# DATA for ZInv
foreach sample( \
DoubleMuParked \
SingleMu \
Photon \
Data_SinglePhoton_Run2012C \
Data_SinglePhoton_Run2012B \
Data_SinglePhotonParked_Run2012D \
) 
  set inputfilelist="/afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set filesPerJob = 1
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob

  foreach jobnumber(`seq 0 1 $maxjob`)
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/RazorControlRegions_${jobnumber}.out -J RazorAnalysis_RazorPhoton_${jobnumber} /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/scripts/runRazorJob_CERN.csh razorPhotonStudy $inputfilelist 1 $filesPerJob $jobnumber RazorControlRegions_${sample}.Job${jobnumber}Of${maxjob}.root /afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/jobs/
    sleep 0.1
  end

end


## CRAB submit    
foreach sample( \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_ext \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph \
ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_ext \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph \
DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext \
WJetsToLNu_HT-150To200_8TeV-madgraph \
WJetsToLNu_HT-200To250_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph \
WJetsToLNu_HT-250To300_8TeV-madgraph_v2 \
WJetsToLNu_HT-300To400_8TeV-madgraph \
WJetsToLNu_HT-300To400_8TeV-madgraph_v2 \
WJetsToLNu_HT-400ToInf_8TeV-madgraph \
WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2 \
QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_20_30_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_80_170_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_170_250_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_250_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
QCD_Pt_350_EMEnriched_TuneZ2star_8TeV_pythia6 \
GJets_HT-40To100_8TeV-madgraph \
GJets_HT-100To200_8TeV-madgraph \
GJets_HT-200To400_8TeV-madgraph_v2 \
GJets_HT-400ToInf_8TeV-madgraph_v3 \
TTGJets_8TeV-madgraph \	
DoubleMuParked \
SingleMu \
Photon \
Data_SinglePhoton_Run2012B \
Data_SinglePhoton_Run2012C \
Data_SinglePhotonParked_Run2012D \
)
  echo "Sample " $sample
  set inputfilelist="/afs/cern.ch/work/a/apresyan/CMSSW_5_3_26/src/RazorAnalyzer/lists/razorNtuplerV1p6-Run1/${sample}.cern.txt"
  set njobs = `cat $inputfilelist | wc | awk '{print $1}' `
  set anatype = razorPhotonStudy
  set option = 0

  sed "s/sampleName/$sample/" crab_runRazorRun.py > crab_tmp.py
  sed -i "s/runRazorCrab/tmp_runRazorCrab/" crab_tmp.py
  sed -i "s/ntupleName.root/RazorAnalysis_$sample.root/" crab_tmp.py 
  sed -i "s|listfile.txt|$inputfilelist|" crab_tmp.py 
  sed -i "s/999/$njobs/" crab_tmp.py
  sed -i "s/999/$njobs/" crab_tmp.py

  sed "s/listfile.txt/${sample}.cern.txt/" runRazorCrab.sh > tmp_runRazorCrab.sh
  sed -i "s/ntupleName.root/RazorAnalysis_$sample.root/" tmp_runRazorCrab.sh
  sed -i "s/Analyses/$anatype/" tmp_runRazorCrab.sh
  sed -i "s/Option/$option/" tmp_runRazorCrab.sh

  crab submit -c crab_tmp.py
  rm crab_tmp.py tmp_runRazorCrab.sh

end

    

foreach sample ( \
Data_SinglePhotonParked_Run2012D \
)
    set inputfilelist = tmp/cms/store/group/phys_susy/razor/Run2Analysis/RazorNtupleV1.6/Run1/Test/MinBias/crab_$sample/
    set stringList = ""
    foreach dir(`find $inputfilelist -type f -name '*root*' | sed -r 's|/[^/]+$||' | sort | uniq`)
    set root = "/root"
    set stringList = "$stringList $dir$root"
    end

echo $stringList | sed  "s/\/root/\/*root/g" > tmp.txt 
set files = `cat tmp.txt`
hadd $sample.root $files && rm tmp.txt
end


