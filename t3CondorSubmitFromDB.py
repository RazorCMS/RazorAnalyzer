#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  t3CondorSubmitRazorAnalzyer.py list1.txt ...
#######################################

analysis = 'razor'
isData = False;
scramArch = "slc5_amd64_gcc462"
cmsswVersion = "CMSSW_5_3_9"

#get the location of this script
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname)
#add RazorProduction to the path and import task2files
razorProductionPath = fullpath+"/../RazorProduction"
if not razorProductionPath in sys.path:
    sys.path.insert(1, razorProductionPath)
import task2files

process = "condor"
datasets = [sys.argv[i] for i in range(1, len(sys.argv))]

#initialize connection to the DB
t2f = task2files.task2files(localdir="/store/group/phys_susy/razor/")

submitScript = open("submitTheCondorJobs.condor", 'w')
submitScript.write('Universe = vanilla\n')
submitScript.write('requirements = Name != "slot1@t3-higgs.ultralight.org" && Name != "slot2@t3-higgs.ultralight.org" && Name != "slot3@t3-higgs.ultralight.org" && Name != "slot4@t3-higgs.ultralight.org" && Name != "slot5@t3-higgs.ultralight.org" && Name != "slot6@t3-higgs.ultralight.org" && Name != "slot7@t3-higgs.ultralight.org" && Name != "slot8@t3-higgs.ultralight.org"\n')
submitScript.write('getenv = True\n\n')

filesperjob = 5
for dataset in datasets:
    output = dataset

    #count lines in file
    with open(inputlist) as f:
        numfiles = sum(1 for _ in f)
    if numfiles == 0: continue
    ijobmax = numfiles/filesperjob + (numfiles % filesperjob > 0)
    extrafiles  = numfiles%ijobmax

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
        outputname = process+"/"+output+"/src/submit_"+str(ijob)+".sh"
        print "OUTPUTFILE: ", outputname
        basedir = fullpath+"/"+process+"/"+output+"/log/";
        outputfile = open(outputname,'w')
        outputfile.write('#!/bin/sh\n')
        outputfile.write('export SCRAM_ARCH='+scramArch+'\n')
        outputfile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        outputfile.write('cd /home/djanders/'+cmsswVersion+'/src\n')
        outputfile.write('eval `scramv1 runtime -sh`\n')
        outputfile.write('cd - 1>/dev/null 2>/dev/null\n')
        outputfile.write('echo $HOSTNAME\n')
        outputfile.write('echo "Current PWD = `pwd`"\n')
        outputfile.write('cp -r '+fullpath+' $PWD/\n')
        outputfile.write('cd $PWD/RazorAnalyzer\n')
        outputfile.write('make clean\n')
        outputfile.write('make\n')
        outputfile.write('./RazorRun '+inputfilename+' '+analysis+' '+fullpath+'/'+process+'/'+output+'/out/'+output+'_'+str(ijob)+'.root\n')
        outputfile.close()
    #    Condor job
        submitScript.write('\nExecutable = '+fullpath+'/'+outputname+'\n')
        submitScript.write('Output = '+process+"/"+output+"/log/"+output+str(ijob)+'.out.$(Process)\n')
        submitScript.write('Log = '+process+"/"+output+"/log/"+str(ijob)+'.log\n')
        submitScript.write('Error = '+process+"/"+output+"/log/"+str(ijob)+'.err\n')
        submitScript.write('Queue\n')
