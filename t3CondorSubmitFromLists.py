#! /usr/bin/env python
import os
import sys, getopt
# set parameters to use cmst3 batch 
#######################################
### usage  t3CondorSubmitRazorAnalzyer.py list1.txt ...
#######################################

isData = False;
scramArch = "slc5_amd64_gcc481"
cmsswVersion = "CMSSW_7_3_0_pre1"
user = os.environ['USER']

#get the location of this script
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname)

if len(sys.argv) < 3:
    print("Usage: t3CondorSubmitFromLists.py analysisName inputList1 inputList2 ...")
    sys.exit()

analysis = sys.argv[1]

process = "condor"
inputlists =  [sys.argv[i] for i in range(2, len(sys.argv))]
datasets = [sys.argv[i].split("/")[-1].replace(".txt","") for i in range(2, len(sys.argv))]

submitScript = open("submitTheCondorJobs.condor", 'w')
submitScript.write('Universe = vanilla\n')
submitScript.write('requirements = Name != "slot1@t3-higgs.ultralight.org" && Name != "slot2@t3-higgs.ultralight.org" && Name != "slot3@t3-higgs.ultralight.org" && Name != "slot4@t3-higgs.ultralight.org" && Name != "slot5@t3-higgs.ultralight.org" && Name != "slot6@t3-higgs.ultralight.org" && Name != "slot7@t3-higgs.ultralight.org" && Name != "slot8@t3-higgs.ultralight.org"\n')
submitScript.write('getenv = True\n\n')

filesperjob = 5
for index, inputlist in enumerate(inputlists):
    output = datasets[index]

    #count lines in file
    with open(inputlist) as f:
        numfiles = sum(1 for _ in f)
    if numfiles == 0: continue
    #ijobmax = 1 #each dataset is a single job 
    ijobmax = numfiles/filesperjob + (numfiles % filesperjob > 0)

    #create directories
    ################################################
    os.system("mkdir -p "+process+"/"+output)
    os.system("mkdir -p "+process+"/"+output+"/log/")
    os.system("mkdir -p "+process+"/"+output+"/input/")
    os.system("mkdir -p "+process+"/"+output+"/src/")
    os.system("mkdir -p "+process+"/"+output+"/out/")
    #######################################
    input = open(inputlist)
    ######################################

    for ijob in range(ijobmax):
        # prepare the list file
        inputfilename = fullpath+"/"+process+"/"+output+"/input/input_"+str(ijob)+".list"
        inputfile = open(inputfilename,'w')
        # if it is a normal job get filesperjob lines
        if ijob != (ijobmax-1):
            for line in range(filesperjob):
                ntpfile = input.readline() 
                inputfile.write(ntpfile)
        else:
            # if it is the last job get ALL remaining lines
            ntpfile = input.readline()
            while ntpfile != '':
                inputfile.write(ntpfile)
                ntpfile = input.readline()
        inputfile.close()

        # prepare the script to run
        outputname = process+"/"+output+"/src/submit_"+output+str(ijob)+".sh"
        print outputname
        basedir = fullpath+"/"+process+"/"+output+"/log/";
        outputfile = open(outputname,'w')
        outputfile.write('#!/bin/sh\n')
        outputfile.write('export JOBDIR=$PWD\n')
        outputfile.write('echo $JOBDIR\n')
        outputfile.write('export SCRAM_ARCH='+scramArch+'\n')
        outputfile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        #outputfile.write('cmsrel '+cmsswVersion+'\n')
        #outputfile.write('cd '+cmsswVersion+'/src\n')
        #outputfile.write('cmsenv\n')
        outputfile.write('cd '+fullpath+'\n')
        outputfile.write('cmsenv\n')
        outputfile.write('cd - \n')
        outputfile.write('echo $HOSTNAME\n')
        outputfile.write(fullpath+'/RazorRun '+inputfilename+' '+analysis+' $JOBDIR/'+analysis+output+str(ijob)+'.root\n')
        #check whether we are on t3-higgs
        outputfile.write('ls /mnt/hadoop/store/user/'+user+' 1>/dev/null 2>/dev/null\n')
        outputfile.write('if [ $? -eq 0 ]\nthen\n   echo "Copying file to hadoop under /store/user/'+user+'/'+process+'/'+output+'"\n   mkdir -p /mnt/hadoop/store/user/'+user+'/'+process+'/'+output+'\n   cp $JOBDIR/'+analysis+output+str(ijob)+'.root /mnt/hadoop/store/user/'+user+'/'+process+'/'+output+'\nelse\n   echo "Copying file to the working directory"\n   cp $JOBDIR/'+analysis+output+str(ijob)+'.root '+fullpath+'/'+process+'/'+output+'/out/\nfi\n')
        outputfile.close()
    #    Condor job
        submitScript.write('\nExecutable = '+fullpath+'/'+outputname+'\n')
        submitScript.write('Output = '+process+"/"+output+"/log/"+output+str(ijob)+'.out\n')
        submitScript.write('Log = '+process+"/"+output+"/log/"+output+str(ijob)+'.log\n')
        submitScript.write('Error = '+process+"/"+output+"/log/"+output+str(ijob)+'.err\n')
        submitScript.write('Queue\n')
