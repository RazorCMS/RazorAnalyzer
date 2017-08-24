#!/bin/sh

# This run script wraps RazorRun_NoAFS and is 
# designed to be executed via HTCondor. It 
# does not read any files from AFS and interacts
# with EOS only to access the input files and
# stage out the output file.

hostname

echo "Arguments: $*"
analysisType=$1
inputfilelist=$2
isData=$3
option=$4
filePerJob=$5
jobnumber=$6
outputfile=$7
outputDirectory=$8
cmsswVersion=$9
label=${10}

echo " "; echo "Initialize CMSSW"; echo " "
workDir=$(pwd)

export SCRAM_ARCH slc6_amd64_gcc530
scramv1 project CMSSW $cmsswVersion
cd $cmsswVersion/src
eval $(scramv1 runtime -sh)
cd -

pwd

echo " "; echo "Show where we are"; echo " "
hostname
pwd

klist

#Do Job splitting and make input file list
cat $inputfilelist | awk "NR > ($jobnumber*$filePerJob) && NR <= (($jobnumber+1)*$filePerJob)" > inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""

datastring=""
if [ ${isData} == 1 ]
then
    datastring="--isData "
fi

# Get ready to run in your home directory
echo " "; echo "Starting razor run job now"; echo " ";
echo ./RazorRun_NoAFS inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring}-f=${outputfile} -n=${option} -l=${label}
./RazorRun_NoAFS inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring}-f=${outputfile} -n=${option} -l=${label} 

ls -ltr 

echo $outputfile 
echo $outputDirectory

#Do below only for output to CERN EOS
eos mkdir -p $outputDirectory
eos cp $outputfile /eos/cms/$outputDirectory/
echo "/eos/cms/$outputDirectory/$outputfile"

status=`echo $?`
echo "Status: $status"

hostname

exit $status
