#!/bin/sh

hostname
echo "Job started"
date

code_dir_suffix=$1
outputDirectory=$2
analysisType=$3
inputfilelist=$4
isData=$5
option=$6
filePerJob=$7
jobnumber=$8
outputfile=$9

currentDir=`pwd`
homeDir=/data/zhicaiz/
runDir=${currentDir}/zhicaiz_${code_dir_suffix}/
rm -rf ${runDir}
mkdir -p ${runDir}

#setup cmssw
cd ${runDir}
echo "entering directory: ${runDir}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_9_4_9
cd CMSSW_9_4_9/src/
export SCRAM_ARCH=slc7_amd64_gcc630
ulimit -c 0
eval `scram runtime -sh`
echo `which root`

git clone https://github.com/RazorCMS/RazorAnalyzer.git
cd RazorAnalyzer/

#get grid proxy
export X509_USER_PROXY=${homeDir}x509_proxy

#run the job
cat ${CMSSW_BASE}${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""

#copy file to local directory
nRemote=`cat inputfilelistForThisJob_${jobnumber}.txt | wc -l`
echo "copy file to local directory ======"

for ifile in `cat inputfilelistForThisJob_${jobnumber}.txt`
do
	xrdcp ${ifile} ./
done

ls razorNtuple_*.root > inputfilelistForThisJob_${jobnumber}_local.txt
nLocal=`cat inputfilelistForThisJob_${jobnumber}_local.txt |wc -l`

if [ "${nRemote}" -eq "${nLocal}" ]
then

	echo "working in directory: `pwd` "
	make
	echo "RazorAnalyzer compiled"
	date
	echo " "; echo "Starting razor run job now"; echo " ";
	echo ./RazorRun_T2 inputfilelistForThisJob_${jobnumber}_local.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile}

	./RazorRun_T2 inputfilelistForThisJob_${jobnumber}_local.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile}

	echo ${outputfile}
	echo ${outputDirectory}
	mkdir -p /mnt/hadoop/${outputDirectory}

	##^_^##
	echo "RazorRun_T2 finished"
	date

	sleep 2
	echo "I slept for 2 second" 

	##job finished, copy file to T2
	echo "copying output file to /mnt/hadoop/${outputDirectory}"
	cp ${outputfile} /mnt/hadoop/${outputDirectory}
	if [ -f /mnt/hadoop/${outputDirectory}/${outputfile} ]
	then
		echo "ZZZZAAAA ============ good news, job finished successfully "
	else
		echo "YYYYZZZZ ============ somehow job failed, please consider resubmitting"
	fi

else
	echo "XXXXYYYY ============= copy file failed, job abandoned"
fi


cd ${currentDir}
rm -rf ${runDir}
echo "Job finished"
date
