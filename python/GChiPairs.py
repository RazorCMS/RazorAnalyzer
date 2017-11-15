#! /usr/bin/env python
import os
import glob

razorSignalDirs = {
        "Razor2016_MoriondRereco": "/eos/cms/store/group/phys_susy/razor/Run2Analysis/FullRazorInclusive/2016/V3p15_05Oct2017/SignalFastsim/"
        }

def parsePair(f):
    pair = f.replace('.root','').split('_')[-2:]
    result = int(pair[0]), int(pair[1])
    return result

def getGChiPairs(signalDir, model):
    pattern = "%s/SMS-%s_*_*.root"%(signalDir, model)
    files = glob.glob(pattern)
    pairs = []
    for f in files:
        # Some model names contain other model names.
        # We check that the filename is actually from
        # the model of interest, 
        mass1Str = os.path.basename(f).replace(
                'SMS-{}_'.format(model), '').split('_')[0]
        try:
            int(mass1Str)
        except ValueError:
            continue
        pairs.append(parsePair(f))
    return pairs

def gchipairs(model, tag='Razor2016_MoriondRereco'):
    gchilist = []
    
    if 'T1x' in model:
        model = 'T1ttbb'
    signalDir = razorSignalDirs[tag].replace('root://eoscms://','')
    return getGChiPairs(signalDir, model)
