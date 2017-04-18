#! /usr/bin/env python
import glob

from macro.razorAnalysis import razorSignalDirs

def parsePair(f):
    pair = f.replace('.root','').split('_')[-2:]
    return int(pair[0]), int(pair[1])

def getGChiPairs(signalDir, model):
    pattern = "%s/%s/SMS-*_*_*.root"%(signalDir, model)
    files = glob.glob(pattern)
    pairs = []
    for f in files:
        pairs.append(parsePair(f))
    return pairs

def gchipairs(model, tag='Razor2016_MoriondRereco'):
    gchilist = []
    
    if 'T1x' in model:
        model = 'T1ttbb'
    signalDir = razorSignalDirs[tag].replace('root://eoscms://','')
    return getGChiPairs(signalDir, model)
