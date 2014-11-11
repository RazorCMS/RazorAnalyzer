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

#get the location of this script
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname)

process = "condor"
inputlists =  [sys.argv[i] for i in range(1, len(sys.argv))]
datasets = [sys.argv[i].split("/")[-1].replace(".txt","") for i in range(1, len(sys.argv))]

submitScript = open("submitTheCondorJobs.sh", 'w')
submitScript.write('#!/bin/sh\n')
filesperjob = 5
for index, inputlist in enumerate(inputlists):
    output = datasets[index]

    #count lines in file
    with open(inputlist) as f:
        numfiles = sum(1 for _ in f)
    ijobmax = numfiles/filesperjob + (numfiles % filesperjob > 0)
    extrafiles  = numfiles%ijobmax

    #create directories
    ################################################
    os.system("mkdir -p "+process+"/"+output)
    os.system("mkdir -p "+process+"/"+output+"/log/")
    os.system("mkdir -p "+process+"/"+output+"/input/")
    os.system("mkdir -p "+process+"/"+output+"/src/")
    os.system("mkdir -p "+process+"/"+output+"/out/")
    os.system("mkdir -p "+process+"/"+output+"/condor/")
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
        outputfile.write('export PATH=/cms/sw/bin:/cms/swslc5/bin:${PATH}\n')
        outputfile.write('export CMS_PATH=/cms/sw\n')
        outputfile.write('cd '+fullpath+'\n')
        outputfile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        outputfile.write('cmsenv\n')
        outputfile.write('make clean\n')
        outputfile.write('make\n')
        outputfile.write('./RazorRun '+inputfilename+' '+analysis+' '+process+'/'+output+'/out/'+output+'_'+str(ijob)+'.root \n')
        outputfile.close()
    #    Condor modification
        submitScript.write("sleep 0.5; condor_submit "+process+"/"+output+"/condor/condor_"+str(ijob)+".sh\n")
        condorScript = open(process+"/"+output+"/condor/condor_"+str(ijob)+".sh", 'w')
        condorScript.write('Executable = '+fullpath+'/'+outputname+'\n')
        condorScript.write('Universe = vanilla\n')
        condorScript.write('Output = '+process+"/"+output+"/out/"+output+str(ijob)+'.out.$(Process)\n')
        condorScript.write('Log = '+process+"/"+output+"/log/"+str(ijob)+'.log\n')
        condorScript.write('Error = '+process+"/"+output+"/log/"+str(ijob)+'.err\n')
        condorScript.write('getenv = True\n\n')
        condorScript.write('Queue\n')
        condorScript.close()
