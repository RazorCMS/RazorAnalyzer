import sys
import os
import argparse
import ROOT as rt

from RunCombine import exec_me

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("inDir", help="input path")
    parser.add_argument("-s", "--smsName", default="T1bbbb", help="SMS name, e.g. T1bbbb")
    parser.add_argument("--dryRun", action="store_true")
    args = parser.parse_args()

    inDir=args.inDir
    smsName=args.smsName
    dryRun=args.dryRun

    #get list of input files
    inFiles = os.listdir(inDir)

    #build dict of files associated with the different signal mass points
    fileLists = {}
    for f in inFiles:
        #parse filename to get gluino and LSP masses
        splitF = f.replace('.root','').split('_')
        mGluino = splitF[-2]
        mLSP = splitF[-1]
        pair = (mGluino, mLSP)

        #add to dictionary if not present
        if pair not in fileLists:
            fileLists[pair] = []

        #add this file to appropriate list
        fileLists[pair].append(f)
        print "Adding",f,"to list of files for model point",pair

    #output directory
    outDir = inDir+'/combined'
    exec_me('mkdir '+outDir, dryRun)
    
    #hadd the files for each signal mass point
    for pair in fileLists:
        print "Signal:",smsName,pair
        exec_me('hadd -f '+outDir+'/SMS-'+smsName+'_'+pair[0]+'_'+pair[1]+'.root '+' '.join([inDir+'/'+f+' ' for f in fileLists[pair]]), dryRun)
