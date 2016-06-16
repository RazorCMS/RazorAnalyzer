#!/bin/env python
import math
import os
import glob
import argparse
from subprocess import call

from ControlRegionNtuples2016 import SAMPLES, SKIMS, DIRS, OPTIONS, VERSION

def submitJobs(analyzer,tag):
    # hard-coded parameters
    queue = '8nh'
    basedir = os.environ['CMSSW_BASE']+'/src/RazorAnalyzer'
    listdir = 'lists/Run2/razorNtupler'+VERSION+'/MC'
    filesperjob = 10
    script=basedir+'/scripts/runRazorJob_CERN_EOS_Dustin.csh'
    os.environ['LSB_JOB_REPORT_MAIL'] = 'N'
    #samples loop
    call(['mkdir','-p',DIRS[tag]+'/jobs'])
    for process in SAMPLES[tag]:
        for sample in SAMPLES[tag][process]:
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
                if not os.path.isfile( os.path.join(DIRS[tag],outfile) ):
                    print "Job %d of %d"%(ijob,maxjob)
                    logfile = os.path.join(basedir,'output','%s_%s_%d.out'%(analyzer,sample,ijob))
                    jobname = '_'.join([analyzer,sample,str(ijob)])
                    cmd = ['bsub','-q',queue,'-o',logfile,'-J',jobname,script,analyzer,inlist,str(0),
                        str(OPTIONS[tag]),str(filesperjob),str(ijob),outfile,
                        DIRS[tag].replace('eos/cms','')+'/jobs', os.environ['CMSSW_BASE']+'/src']
                    print ' '.join(cmd)
                    call(cmd)

def haddFiles(analyzer,tag):
    for process in SAMPLES[tag]:
        for sample in SAMPLES[tag][process]:
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tag',help='1L, 2L, ...')
    parser.add_argument('--analyzer',default='RazorControlRegions',help='Analyzer name')
    parser.add_argument('--submit',action='store_true',help='Submit batch jobs')
    parser.add_argument('--hadd',action='store_true',help='Combine ntuple files')
    parser.add_argument('--normalize',action='store_true',help='Normalize ntuple files')
    parser.add_argument('--hadd-final',dest='haddFinal',action='store_true',help='Combine normalized ntuple files')
    args = parser.parse_args()
    analyzer = args.analyzer
    tag = args.tag

    print "Analyzer:",analyzer
    print "Tag:",tag

    if args.submit:
        print "Submit batch jobs..."
        submitJobs(analyzer,tag)

    if args.hadd:
        print "Combine ntuples..."
        haddFiles(analyzer,tag)

    if args.normalize:
        print "Normalize ntuples..."
        normalizeFiles(analyzer,tag)

    if args.haddFinal:
        print "Combine normalized files..."
        haddNormalizedFiles(analyzer,tag)
