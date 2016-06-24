#!/bin/env python
import math
import os,sys
import glob
import argparse
from subprocess import call, check_output

from ControlRegionNtuples2016 import SAMPLES, SKIMS, DIRS, OPTIONS, VERSION, DATA

def submitJobs(analyzer,tag,isData=False,submit=False):
    # parameters
    samples = SAMPLES
    queue = '8nh'
    basedir = os.environ['CMSSW_BASE']+'/src/RazorAnalyzer'
    listdir = 'lists/Run2/razorNtupler'+VERSION+'/MC'
    if isData:
        listdir = listdir.replace('/MC','/data')
        samples = DATA
    filesperjob = 10
    script=basedir+'/scripts/runRazorJob_CERN_EOS_Dustin.csh'
    os.environ['LSB_JOB_REPORT_MAIL'] = 'N'
    #samples loop
    call(['mkdir','-p',DIRS[tag]+'/jobs'])
    for process in samples[tag]:
        for sample in samples[tag][process]:
            inlist = os.path.join(basedir,listdir,sample+'.cern.txt')
            if not os.path.isfile(inlist):
                print "Warning: list file not found for",sample
                continue
            nfiles = sum([1 for line in open(inlist)])
            maxjob = int(math.ceil( nfiles*1.0/filesperjob ))-1
            print "Sample:",sample," maxjob =",maxjob
            #submit
            for ijob in range(maxjob+1):
                outfile = '%s_%s.Job%dof%d.root'%(analyzer,sample,ijob,maxjob)
                if not os.path.isfile( os.path.join(DIRS[tag],'jobs',outfile) ):
                    print "Job %d of %d"%(ijob,maxjob)
                    logfile = os.path.join(basedir,'output','%s_%s_%d.out'%(analyzer,sample,ijob))
                    jobname = '_'.join([analyzer,sample,str(ijob)])
                    cmd = ['bsub','-q',queue,'-o',logfile,'-J',jobname,script,analyzer,inlist,
                            str(int(isData)),str(OPTIONS[tag]),str(filesperjob),str(ijob),outfile,
                        DIRS[tag].replace('eos/cms','')+'/jobs', os.environ['CMSSW_BASE']+'/src']
                    print ' '.join(cmd)
                    if submit:
                        call(cmd)

def haddFiles(analyzer,tag,isData=False):
    samples = SAMPLES
    if isData:
        samples = DATA
    for process in samples[tag]:
        for sample in samples[tag][process]:
            print "Sample:",sample
            fname = DIRS[tag]+'/'+analyzer+'_'+sample+'.root'
            query = DIRS[tag]+'/jobs/'+analyzer+'_'+sample+'.Job*.root'
            jobfiles = glob.glob( query )
            if os.path.isfile( fname ):
                print "File",fname,"exists; skipping"
            elif len(jobfiles) > 0:
                call(['hadd',fname]+jobfiles)
            else:
                print "Warning: no files found (",query,")"

def normalizeFiles(analyzer,tag):
    #make list file for normalizing
    paths = glob.glob( DIRS[tag]+'/*.root' )
    with open('ntuples_'+tag+'.txt','w') as normfile:
        for f in paths:
            sample = os.path.basename(f).replace('.root','').replace(analyzer+'_','')
            #check if we need this sample
            for process in SAMPLES[tag]:
                if sample in SAMPLES[tag][process]:
                    normfile.write(sample+' '+f+'\n')
                    break
    call(['./NormalizeNtuple','ntuples_'+tag+'.txt'])

def haddNormalizedFiles(analyzer,tag):
    for process in SAMPLES[tag]:
        print "Process:",process
        haddList = []
        for sample in SAMPLES[tag][process]:
            print "Sample:",sample,
            fname = DIRS[tag]+'/'+analyzer+'_'+sample+'_1pb_weighted.root'
            if os.path.isfile( fname ):
                print "found!"
                haddList.append( fname )
            else:
                print "not found..."
        outName = DIRS[tag]+'/'+analyzer+'_'+process+'_1pb_weighted.root'
        if len(haddList) > 0:
            call(['hadd',outName]+haddList)
        else:
            print "Skipping process",process,"because no input files were found!"

def removeDuplicates(analyzer,tag):
    #first hadd all datasets together
    haddList = []
    for process in DATA[tag]:
        print "Dataset:",process
        for sample in DATA[tag][process]:
            fname = DIRS[tag]+'/'+analyzer+'_'+sample+'.root'
            print fname,
            if os.path.isfile( fname ):
                print "found!"
                haddList.append( fname )
            else:
                print "not found..."
    outName = DIRS[tag]+'/'+analyzer+'_Data.root'
    if len(haddList) > 0:
        call(['hadd',outName]+haddList)
    else:
        print "No input files were found!"
        return
    #now remove duplicates
    outNameNoDuplicates = outName.replace('_Data','_Data_NoDuplicates')
    dupMacro = 'macros/RemoveDuplicateEvents.C'
    macroCall = "root -l '%s++(\"%s\",\"%s\")'"%(dupMacro,outName,outNameNoDuplicates)
    print macroCall
    call(macroCall,shell=True)

def goodLumi(analyzer,tag):
    #get location of good lumi script and python file containing JSON link
    script = check_output(['which', 'FWLiteGoodLumi'])
    pythonFile = os.path.join(os.environ['CMSSW_BASE'],'src','RazorCommon','Tools','python','loadJson.py')
    #copy file locally
    localInName = analyzer+'_Data_NoDuplicates.root'
    localOutName = localInName.replace('.root','_GoodLumiGolden.root')
    inName = DIRS[tag]+'/'+localInName
    call(['cp',inName,'.'])
    #call script
    print "FWLiteGoodLumi %s %s %s"%(pythonFile,localInName,localOutName)
    call(['FWLiteGoodLumi',pythonFile,localInName,localOutName])
    #move to target directory
    call(['mv',localOutName,DIRS[tag]])
    #remove temp file
    call(['rm',localInName])

def copyLocal(analyzer,tag,isData=False):
    samples = SAMPLES
    if isData:
        samples = DATA
    #make directory
    localdir = 'Backgrounds/'+tag
    call(['mkdir','-p',localdir])
    if isData:
        fname = DIRS[tag]+'/'+analyzer+'_Data_NoDuplicates_GoodLumiGolden.root'
        print "cp",fname,localdir
        call(['cp',fname,localdir])
    else:
        for process in samples[tag]:
            fname = DIRS[tag]+'/'+analyzer+'_'+process+'_1pb_weighted.root'
            print "cp",fname,localdir
            call(['cp',fname,localdir])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag',help='1L, 2L, ...')
    parser.add_argument('--analyzer',default='RazorControlRegions',help='Analyzer name')
    parser.add_argument('--data',action='store_true',help='Run on data (MC otherwise)')
    parser.add_argument('--submit',action='store_true',help='Submit batch jobs')
    parser.add_argument('--no-sub', dest='noSub',action='store_true', 
            help='Print commands but do not submit')
    parser.add_argument('--hadd',action='store_true',help='Combine ntuple files')
    parser.add_argument('--normalize',action='store_true',help='Normalize ntuple files')
    parser.add_argument('--hadd-final',dest='haddFinal',action='store_true',help='Combine normalized ntuple files')
    parser.add_argument('--remove-duplicates',dest='removeDuplicates',action='store_true',help='Remove duplicates')
    parser.add_argument('--good-lumi',dest='goodLumi',action='store_true',help='Apply good lumi selection')
    parser.add_argument('--copy-local',dest='copyLocal',action='store_true',help='Copy files locally')
    args = parser.parse_args()
    analyzer = args.analyzer
    tag = args.tag
    isData = args.data
    noSub = args.noSub

    #check if EOS is mounted
    if not os.path.isdir('eos/cms/store'):
        sys.exit("Please mount EOS under ./eos before using this tool.")

    print "Analyzer:",analyzer
    print "Tag:",tag

    if args.submit:
        print "Submit batch jobs..."
        submitJobs(analyzer,tag,isData,submit=(not noSub))

    if args.hadd:
        print "Combine ntuples..."
        haddFiles(analyzer,tag,isData)

    if args.normalize:
        print "Normalize ntuples..."
        if isData:
            sys.exit("Error: options --data and --normalize do not make sense together!")
        normalizeFiles(analyzer,tag)

    if args.haddFinal:
        print "Combine normalized files..."
        if isData:
            sys.exit("Error: options --data and --hadd-final do not make sense together!")
        haddNormalizedFiles(analyzer,tag)

    if args.removeDuplicates:
        print "Remove duplicate events..."
        if not isData:
            print "--data was not specified, but I assume you want to use data."
        removeDuplicates(analyzer,tag)

    if args.goodLumi:
        print "Apply good lumi selection..."
        if not isData:
            print "--data was not specified, but I assume you want to use data."
        goodLumi(analyzer,tag)

    if args.copyLocal:
        print "Copy files locally..."
        copyLocal(analyzer,tag,isData)
