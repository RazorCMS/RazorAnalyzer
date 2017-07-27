#!/bin/tcsh

hostname

echo "Arguments: $*"
set analysisType=$1
set inputfilelist=$2
set isData=$3
set option=$4
set filePerJob=$5
set jobnumber=$6
set outputfile=$7
set outputDirectory=$8
set cmsswDir=$9
set label=$10

echo " "; echo "Initialize CMSSW"; echo " "
set workDir=`pwd`

setenv SCRAM_ARCH slc6_amd64_gcc481
cd    $cmsswDir
eval `scramv1 runtime -csh`
cd -

pwd

cp $CMSSW_BASE/src/RazorAnalyzer/RazorRun ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSVv2_Moriond17_B_H.csv ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSVv2_Moriond17_G_H.csv ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSVv2_ichep.csv ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/fastsim_csvv2_ttbar_26_1_2017.csv ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSV_13TEV_Combined_20_11_2015.csv ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/JEC/JEC_Summer16_23Sep2016V3.tgz ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/JEC/Spring16_FastSimV1.tgz ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/JEC/JetResolutionInputAK5PF.txt ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonCorrections/Winter_2016_reReco_v1_ele_smearings.dat ./
cp -v /eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonCorrections/Winter_2016_reReco_v1_ele_scales.dat ./
tar vxzf JEC_Summer16_23Sep2016V3.tgz
tar vxzf Spring16_FastSimV1.tgz


echo " "; echo "Show where we are"; echo " "
hostname
pwd

klist

#Do Job splitting and make input file list
cat $inputfilelist | awk "NR > ($jobnumber*$filePerJob) && NR <= (($jobnumber+1)*$filePerJob)" >! inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""

set datastring=""
if (${isData} == 1) then
    set datastring="--isData "
endif

# Get ready to run in your home directory
echo " "; echo "Starting razor run job now"; echo " ";
echo ./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring}-f=${outputfile} -n=${option}
./RazorRun inputfilelistForThisJob_${jobnumber}.txt ${analysisType} ${datastring}-f=${outputfile} -n=${option} -l=${label} |& tee ${outputfile}.log

ls -ltr 

echo $outputfile 
echo $outputDirectory

#Do below only for output to CERN EOS
#cmsMkdir $outputDirectory
#cmsStage -f $outputfile $outputDirectory

mkdir -p /eos/cms/$outputDirectory
cp -v $outputfile /eos/cms/$outputDirectory

set status=`echo $?`
echo "Status: $status"

hostname

exit $status
