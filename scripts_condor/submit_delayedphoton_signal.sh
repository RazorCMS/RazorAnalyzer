#!/bin/sh

mkdir -p log
mkdir -p submit

cd ../
RazorAnalyzerDir=`pwd`
cd -

job_script=${RazorAnalyzerDir}/scripts_condor/runRazorJob_CaltechT2.sh
filesPerJob=20

for sample in \
GMSB_L100TeV_Ctau1000cm_13TeV-pythia8 \
GMSB_L100TeV_Ctau10cm_13TeV-pythia8 

do
	echo "Sample " ${sample}
	inputfilelist=${RazorAnalyzerDir}/lists/Run2/razorNtuplerV4p1/MC_Summer16_reMINIAOD/${sample}.cern.txt
	nfiles=`cat $inputfilelist | wc | awk '{print $1}' `
	maxjob=`python -c "print int($nfiles.0/$filesPerJob)-1"`
	rm submit/${sample}_Job*.jdl

	for jobnumber in `seq 0 1 ${maxjob}`
	do
		echo "job " ${jobnumber} " out of " ${maxjob}
		jdl_file=submit/${sample}_Job${jobnumber}_Of_${maxjob}.jdl
		echo "Universe = vanilla" > ${jdl_file}
		echo "Executable = ${job_script}" >> ${jdl_file}
		echo "Arguments = DelayedPhoton_GMSB_L350TeV_Ctau200cm_13TeV-pythia8_Job${jobnumber}_Of_${maxjob} /store/user/zhicaiz/Run2Analysis/DelayedPhotonAnalysis/2016/jobs/ DelayedPhotonAnalyzer ${inputfilelist} no 10 ${filesPerJob} ${jobnumber} DelayedPhoton_GMSB_L350TeV_Ctau200cm_13TeV-pythia8_Job${jobnumber}_Of_${maxjob}.root" >> ${jdl_file}
		echo "Log = log/${sample}_Job${jobnumber}_Of_${maxjob}_PC.log" >> ${jdl_file}
		echo "Output = log/${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).out" >> ${jdl_file}
		echo "Error = log/${sample}_Job${jobnumber}_Of_${maxjob}_\$(Cluster).\$(Process).err" >> ${jdl_file}
		echo "should_transfer_files = YES" >> ${jdl_file}
		echo "RequestMemory = 2000" >> ${jdl_file}
		echo "RequestCpus = 1" >> ${jdl_file}
		echo "when_to_transfer_output = ON_EXIT" >> ${jdl_file}
		echo "Queue 1" >> ${jdl_file}
		echo "condor_submit submit/${sample}_Job${jobnumber}_Of_${maxjob}.jdl"
		condor_submit submit/${sample}_Job${jobnumber}_Of_${maxjob}.jdl
	done
done

