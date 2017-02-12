#!/bin/tcsh

hostname

echo "Arguments: $*"
set analysisType=$1
set inputfilelist=$2
set isData=$3
set option=$4
set label=$5
set filePerJob=$6
set jobnumber=$7
set outputfile=$8
set outputDirectory=$9

echo " "; echo "Initialize CMSSW"; echo " "
#setenv KRB5CCNAME /home/sixie/.krb5/ticket
set workDir=`pwd`

setenv SCRAM_ARCH slc6_amd64_gcc491
cd    /afs/cern.ch/work/s/sixie/public/releases/run2/CMSSW_7_4_2/src/
eval `scramv1 runtime -csh`
cd -

pwd
# env

cp $CMSSW_BASE/src/RazorAnalyzer/bin/Run${analysisType} ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSVv2_Moriond17_B_H.csv ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSVv2_Moriond17_G_H.csv ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSVv2_ichep.csv ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/fastsim_csvv2_ttbar_26_1_2017.csv ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/BTagEfficiencies/CSV_13TEV_Combined_20_11_2015.csv ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/JEC/JEC_Summer16_23Sep2016V3.tgz ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/JEC/Spring16_FastSimV1.tgz ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/JEC/JetResolutionInputAK5PF.txt ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonCorrections/Winter_2016_reReco_v1_ele_smearings.dat ./
eos cp root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PhotonCorrections/Winter_2016_reReco_v1_ele_scales.dat ./
tar vxzf JEC_Summer16_23Sep2016V3.tgz
tar vxzf Spring16_FastSimV1.tgz

echo " "; echo "Show where we are"; echo " "
hostname
pwd
## env

klist

#setenv STAGE_SVCCLASS cmsprod

#Do Job splitting and make input file list
cat $inputfilelist | awk "NR > ($jobnumber*$filePerJob) && NR <= (($jobnumber+1)*$filePerJob)" >! inputfilelistForThisJob_${jobnumber}.txt
echo ""
echo "************************************"
echo "Running on these input files:"
cat inputfilelistForThisJob_${jobnumber}.txt
echo "************************************"
echo ""


# Get ready to run in your home directory
if (${isData} == "false") then
   set isDataOption = ""
else 
   set isDataOption = "--isData"
endif

echo " "; echo "Starting razor run job now"; echo " ";
echo ./Run${analysisType} inputfilelistForThisJob_${jobnumber}.txt ${isDataOption} --outputFile=${outputfile} --optionNumber=${option}
./Run${analysisType} inputfilelistForThisJob_${jobnumber}.txt ${isDataOption} --outputFile=${outputfile} --optionNumber=${option} -l=${label} |& tee ${outputfile}.log

ls -ltr 

echo $outputfile 
echo $outputDirectory

#Do below only for output to CERN EOS
cmsMkdir $outputDirectory
cmsStage -f $outputfile $outputDirectory
#cmsStage -f ${outputfile}.log $outputDirectory

set tempOutputfile = `echo $outputfile | sed 's/.root//'`
foreach f ( ${tempOutputfile}_*.root )
   cmsStage -f $f $outputDirectory
end

#mkdir -p eos
#eosmount eos
#cp -v $outputfile eos/cms/$outputDirectory
#eosumount eos

set status=`echo $?`
echo "Status: $status"

hostname

exit $status
