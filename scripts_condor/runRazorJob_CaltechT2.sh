#!/bin/sh

hostname
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


thisDir=/tmp/zhicaiz_${code_dir_suffix}/
rm -rf ${thisDir}
mkdir -p ${thisDir}

homeDir=/data/zhicaiz/
#setup cmssw
cd ${homeDir}release/RazorAnalyzer/CMSSW_9_4_9/src/RazorAnalyzer/
workDir=`pwd`
echo "entering directory: ${workDir}"
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
ulimit -c 0
eval `scram runtime -sh`

echo `which root`
#get grid proxy
#source /data/zhicaiz/crab.sh
export X509_USER_PROXY=${homeDir}x509_proxy

#copy RazorRun_T2 to local directory
cd ${thisDir}
echo "back to directory: ${thisDir}"
cp $CMSSW_BASE/src/RazorAnalyzer/RazorRun_T2 ./

#run the job
cat ${inputfilelist} | awk "NR > (${jobnumber}*${filePerJob}) && NR <= ((${jobnumber}+1)*${filePerJob})" > inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""

echo " "; echo "Starting razor run job now"; echo " ";
echo ./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile}

./RazorRun_T2 inputfilelistForThisJob_${jobnumber}.txt ${analysisType} -d=${isData} -n=${option} -f=${outputfile}

echo ${outputfile}
echo ${outputDirectory}
mkdir -p /mnt/hadoop/${outputDirectory}

##^_^##
sleep 1
echo "I slept for 1 second" 

##job finished, copy file to T2
echo "copying output file to /mnt/hadoop/${outputDirectory}"
cp ${outputfile} /mnt/hadoop/${outputDirectory}

cd ${workDir}/scripts_condor
date
