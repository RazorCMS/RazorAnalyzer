#!/bin/tcsh

foreach sample( \
WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8 \
#WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8 \
#ZGGJetsToLLGG_5f_LO_amcatnloMLM_pythia8 \
#ZGGToLLGG_5f_TuneCUETP8M1_13TeV-amcatnlo-pythia8 \
#ZGGJets_ZToHadOrNu_5f_LO_madgraph_pythia8 \
#ZGGToNuNuGG_5f_TuneCUETP8M1_13TeV-amcatnlo-pythia8 \
WGG_5f_TuneCUETP8M1_13TeV-amcatnlo-pythia8 \
WGGJets_TuneCUETP8M1_13TeV_madgraphMLM_pythia8 \
TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8 \
ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8 \
WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph \
TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8 \
TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8 \
)
  setenv LSB_JOB_REPORT_MAIL N 
  set inputfilelist="/afs/cern.ch/work/c/cpena/public/MultiBoson/CMSSW_9_2_0/src/RazorAnalyzer/lists/Run2/razorNtuplerV3p8/MC_Summer16/${sample}.cern.txt"
  set filesPerJob = 1 
  set nfiles = `cat $inputfilelist | wc | awk '{print $1}' `
  set maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
  echo "Sample " $sample " maxjob = " $maxjob 
  mkdir -p ../submissionbkg
  foreach jobnumber(`seq 0 1 $maxjob`)
    if ( ! -e /afs/cern.ch/user/s/sixie/eos/cms/store/group/phys_susy/razor/SusyEwkHgg/jobs/SusyEwkHgg_${sample}.Job${jobnumber}Of${maxjob}.root ) then
    echo "job " $jobnumber " out of " $maxjob
    bsub -q 8nm -o /afs/cern.ch/work/c/cpena/public/MultiBoson/CMSSW_9_2_0/src/RazorAnalyzer/submissionbkg/SusyEwkHgg_${sample}_${jobnumber}.out -J SusyEwkHgg_${sample}_${jobnumber} /afs/cern.ch/work/c/cpena/public/MultiBoson/CMSSW_9_2_0/src/RazorAnalyzer/scripts/runRazorJob_CERN_EOS.csh SusyEwkHgg $inputfilelist 0 0 $filesPerJob $jobnumber SusyEwkHgg_${sample}.Job${jobnumber}Of${maxjob}.root /store/group/phys_susy/razor/Run2Analysis/SusyEwkHgg/bkg/cp/${sample}_jobs/ /afs/cern.ch/work/c/cpena/public/MultiBoson/CMSSW_9_2_0/src/RazorAnalyzer/ Razor2016_MoriondRereco
    sleep 0.1
  end

end


