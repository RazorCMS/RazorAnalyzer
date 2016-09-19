#!/bin/env python

### Sequence for data: submit --> hadd --> skim --> hadd-final --> remove-duplicates --> good-lumi --> copy-local
### Sequence for MC: submit --> hadd --> normalize --> hadd-final --> copy-local

import math
import os,sys
import glob
import argparse
from subprocess import call, check_output

from ControlRegionNtuples2016_V3p5 import SAMPLES, TREETYPES, TREETYPEEXT, SKIMS, DIRS, OPTIONS, VERSION, DATA, SUFFIXES

def getSamplePrefix(analyzer,tag,reHLT=False):
    return analyzer.replace('RazorControl','RunTwoRazorControl')+(
            (TREETYPEEXT[tag]!='')*('_'+TREETYPEEXT[tag]))+(
            (SKIMS[tag]!='')*('_'+SKIMS[tag]))+(
            (reHLT==True)*('_reHLT'))

def getJobFileName(analyzer,tag,sample,ijob,maxjob,reHLT=False):
    prefix = getSamplePrefix(analyzer,tag,reHLT)
    return '%s_%s.Job%dof%d.root'%(prefix,sample,ijob,maxjob)

def getFileName(analyzer,tag,sample,reHLT=False):
    prefix = getSamplePrefix(analyzer,tag,reHLT)
    return '%s_%s.root'%(prefix,sample)

def submitJobs(analyzer,tag,isData=False,submit=False,reHLT=False):
    # parameters
    samples = SAMPLES
    queue = '1nh'
    basedir = os.environ['CMSSW_BASE']+'/src/RazorAnalyzer'
    listdir = 'lists/Run2/razorNtupler'+(VERSION.split('_')[0])+'/MC'
    if reHLT: listdir += 'reHLT'
    jobssuffix = '/jobs'
    if isData:
        listdir = listdir.replace('/MC','/data')
        samples = DATA
    filesperjob = 10
    script=basedir+'/scripts/runRazorJob_CERN_EOS_Dustin.csh'
    os.environ['LSB_JOB_REPORT_MAIL'] = 'N'
    #samples loop
    call(['mkdir','-p',DIRS[tag]+'/jobs'])
    print samples[tag]
    for process in samples[tag]:
        print process
        for sample in samples[tag][process]:
            print sample
            inlist = os.path.join(basedir,listdir,sample+'.cern.txt')
            if not os.path.isfile(inlist):
                print "Warning: list file",inlist,"not found!"
                continue
            nfiles = sum([1 for line in open(inlist)])
            maxjob = int(math.ceil( nfiles*1.0/filesperjob ))-1
            print "Sample:",sample," maxjob =",maxjob
            #submit
            for ijob in range(maxjob+1):
                outfile = getJobFileName(analyzer,tag,sample,ijob,maxjob,reHLT)
                if not os.path.isfile( DIRS[tag]+jobssuffix+'/'+outfile ):
                    print "Job %d of %d"%(ijob,maxjob)
                    logfile = os.path.join(basedir,'output','%s_%s_%d.out'%(analyzer,sample,ijob))
                    jobname = '_'.join([analyzer,sample,str(ijob)])
                    cmd = ['bsub','-q',queue,'-o',logfile,'-J',jobname,script,analyzer,inlist,
                            str(int(isData)),str(OPTIONS[tag]),str(filesperjob),str(ijob),outfile,
                        DIRS[tag].replace('eos/cms','')+jobssuffix, os.environ['CMSSW_BASE']+'/src']
                    print ' '.join(cmd)
                    if submit:
                        call(cmd)

def haddFiles(analyzer,tag,isData=False,force=False,reHLT=False):
    samples = SAMPLES
    if isData:
        samples = DATA
    for process in samples[tag]:
        for sample in samples[tag][process]:
            print "Sample:",sample
            fname = DIRS[tag]+'/'+getFileName(analyzer,tag,sample,reHLT)
            query = DIRS[tag]+'/jobs/'+(getFileName(analyzer,tag,sample,reHLT).replace('.root','*.Job*.root'))
            jobfiles = glob.glob( query )
            if os.path.isfile( fname ) and not force:
                print "File",fname,"exists; skipping"
            elif len(jobfiles) > 0:
                if force:
                    call(['hadd','-f',fname]+jobfiles)
                else:
                    call(['hadd',fname]+jobfiles)
            else:
                print "Warning: no files found (",query,")"

def normalizeFiles(analyzer,tag,force=False,reHLT=False):
    #make list file for normalizing
    paths = glob.glob( DIRS[tag]+'/*.root' )
    with open('ntuples_'+tag+'.txt','w') as normfile:
        for f in paths:
            #check if normalized file exists
            if (not force) and os.path.isfile( f.replace('.root','_1pb_weighted.root') ): continue
            sample = os.path.basename(f).replace('.root','').replace(
                    getSamplePrefix(analyzer,tag,reHLT)+'_','')
            #check if we need this sample
            for process in SAMPLES[tag]:
                if sample in SAMPLES[tag][process]:
                    normfile.write(sample+' '+f+'\n')
                    break
    call(['./NormalizeNtuple','ntuples_'+tag+'.txt'])

def haddNormalizedFiles(analyzer,tag,force=False,reHLT=False):
    for process in SAMPLES[tag]:
        print "Process:",process
        haddList = []
        for sample in SAMPLES[tag][process]:
            print "Sample:",sample,
            fname = DIRS[tag]+'/'+getFileName(analyzer,tag,sample,reHLT).replace(
                    '.root','_1pb_weighted.root')
            if os.path.isfile( fname ):
                print "found!"
                haddList.append( fname )
            else:
                print "not found..."
        outName = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag,reHLT)+'_'+process+'_1pb_weighted.root'
        if len(haddList) > 0:
            if force:
                call(['hadd','-f',outName]+haddList)
            else:
                call(['hadd',outName]+haddList)
        else:
            print "Skipping process",process,"because no input files were found!"
            print "Failed to create file",outName

def combineData(analyzer,tag,force=False):
    haddList = []
    for process in DATA[tag]:
        print "Dataset:",process
        for sample in DATA[tag][process]:
            fname = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag)+'_'+sample+'_RazorSkim.root'
            print fname,
            if os.path.isfile( fname ):
                print "found!"
                haddList.append( fname )
            else:
                print "not found..."
    outName = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag)+'_Data_RazorSkim.root'
    if len(haddList) > 0:
        if force:
            call(['hadd','-f',outName]+haddList)
        else:
            call(['hadd',outName]+haddList)
    else:
        print "No input files were found!"
        return

def skimNtuples(analyzer,tag,isData=False):
    samples = SAMPLES
    if isData:
        samples = DATA
    #make list file for skimming
    with open('skim_'+tag+'.txt','w') as skimfile:
        for process in samples[tag]:
            for sample in samples[tag][process]:
                fname = DIRS[tag]+'/'+getFileName(analyzer,tag,sample)
                if not isData: fname = fname.replace('.root','_1pb_weighted.root')
                if os.path.isfile( fname ):
                    print "Adding",sample,"to skim file"
                    skimfile.write(fname+'\n')
                else:
                    print "Input file for",sample,"not found!"
                    print "( looking for",fname,")"
    skimString = 'MR%s > 300 && Rsq%s > 0.15'%(SUFFIXES[tag],SUFFIXES[tag])
    print "Skimming with",skimString
    call(['./SkimNtuple','skim_'+tag+'.txt',DIRS[tag],'RazorSkim',skimString])

def removeDuplicates(analyzer,tag):
    outName = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag)+'_Data_RazorSkim.root'
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
    localInName = getSamplePrefix(analyzer,tag)+'_Data_NoDuplicates_RazorSkim.root'
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
        fname = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag)+'_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root'
        print "cp",fname,localdir
        call(['cp',fname,localdir])
    else:
        for process in samples[tag]:
            fname = DIRS[tag]+'/'+getSamplePrefix(analyzer,tag)+'_'+process+'_1pb_weighted.root'
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
    parser.add_argument('--force', action='store_true', help='Force HADD to recreate files')
    parser.add_argument('--hadd',action='store_true',help='Combine ntuple files')
    parser.add_argument('--normalize',action='store_true',help='Normalize ntuple files')
    parser.add_argument('--hadd-final',dest='haddFinal',action='store_true',help='Combine normalized ntuple files')
    parser.add_argument('--skim',action='store_true',help='Apply razor skim to final ntuples')
    parser.add_argument('--remove-duplicates',dest='removeDuplicates',action='store_true',help='Remove duplicates')
    parser.add_argument('--good-lumi',dest='goodLumi',action='store_true',help='Apply good lumi selection')
    parser.add_argument('--copy-local',dest='copyLocal',action='store_true',help='Copy files locally')
    parser.add_argument('--reHLT',action='store_true',help='Process reHLT samples')
    args = parser.parse_args()
    tag = args.tag
    analyzer = args.analyzer
    if tag == 'Signal':
        analyzer = 'FullRazorInclusive'
    elif tag == 'LiteZMu' or tag == 'LiteZEle':
        analyzer = 'RazorLiteZ'
    isData = args.data
    noSub = args.noSub
    force = args.force
    reHLT = args.reHLT

    #check if EOS is mounted
    if not os.path.isdir('eos/cms/store'):
        sys.exit("Please mount EOS under ./eos before using this tool.")

    print "Analyzer:",analyzer
    print "Tag:",tag

    if args.submit:
        print "Submit batch jobs..."
        submitJobs(analyzer,tag,isData,submit=(not noSub),reHLT=reHLT)

    if args.hadd:
        print "Combine ntuples..."
        haddFiles(analyzer,tag,isData,force,reHLT=reHLT)

    if args.normalize:
        print "Normalize ntuples..."
        if isData:
            sys.exit("Error: options --data and --normalize do not make sense together!")
        normalizeFiles(analyzer,tag,force,reHLT=reHLT)

    if args.haddFinal:
        print "Combine normalized files..."
        if isData:
            combineData(analyzer,tag,force)
        else:
            haddNormalizedFiles(analyzer,tag,force,reHLT=reHLT)

    if args.skim:
        print "Skim finished ntuples..."
        skimNtuples(analyzer,tag,isData)

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
