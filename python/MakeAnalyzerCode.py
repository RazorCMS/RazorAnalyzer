#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print "Please specify an analyzer name!"
    sys.exit()

analyzer = sys.argv[1]

inNames = ['include/AnalyzerTemplate.txt','src/RunAnalyzerTemplate.txt']
outNames = ['analyzers/'+analyzer+'.h','src/Run'+analyzer+'.cc']

for i in range(len(inNames)):
    with open(inNames[i]) as inF:
        with open(outNames[i],'w') as outF:
            for line in inF:
                outF.write( line.replace('%ANALYZER%',analyzer) )       
