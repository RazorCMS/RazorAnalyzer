#!/bin/env python

import os, sys, subprocess
import argparse

def sub_sequence(tag, isData=False, submit=False, label=''):
    basedir = os.environ['CMSSW_BASE']+'/src/RazorAnalyzer'
    if not submit: 
        nosub = '--no-sub'
    else:
        nosub = ''
    if isData:
        data = '--data'
    else:
        data = ''
    

    cmd_submit = list(filter(None,['python', 'python/ntupling/NtupleUtils.py', tag, '--submit', nosub, '--label', label, data]))
    print ' '.join(cmd_submit)
    subprocess.call(cmd_submit)
    if submit: 
        # Get the number of running and pending jobs, only proceed if both = 0
        jstatus = subprocess.Popen(["bjobs","-sum","-noheader","-J","FullRazorInclusive*"], stdout=subprocess.PIPE).communicate()[0]
        jlist = [ x for x in jstatus.split(" ") if x != '' ]
        rjobs = int(jlist[0])
        pjobs = int(jlist[4])
        if rjobs == 0 and pjobs == 0:
            cmd_hadd = list(filter(None,['python', 'python/ntupling/NtupleUtils.py', tag, '--hadd', nosub, '--label', label, data]))
            print ' '.join(cmd_hadd)
            subprocess.call(cmd_hadd)
            if not isData:
                cmd_normalize = list(filter(None,['python', 'python/ntupling/NtupleUtils.py', tag, '--normalize', nosub, '--label', label, data]))
                print ' '.join(cmd_normalize)
                subprocess.call(cmd_normalize)
            cmd_hadd_final = list(filter(None,['python', 'python/ntupling/NtupleUtils.py', tag, '--hadd-final', nosub, '--label', label, data]))
            print ' '.join(cmd_hadd_final)
            subprocess.call(cmd_hadd_final)
            if isData:
                cmd_remove_duplicates = list(filter(None,['python', 'python/ntupling/NtupleUtils.py', '--remove-duplicates', nosub, '--label', label, data, tag]))
                print ' '.join(cmd_remove_duplicates)
                subprocess.call(cmd_remove_duplicates)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tag', help = '1L, 2L, ...', required = True)
    parser.add_argument('--label', help = 'Label for RazorRun', default = '')
    parser.add_argument('--data', action = 'store_true', help = 'Run on data (MC otherwise)')
    parser.add_argument('--no-sub', dest = 'noSub', action = 'store_true', help = 'Print commands but do not submit')

    args = parser.parse_args()
    
    sub_sequence(args.tag, args.data, (not args.noSub), args.label)
