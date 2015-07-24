import json
from optparse import OptionParser
import ROOT as rt
import sys
from array import *
import os
from itertools import *
from operator import *
import pickle

def walk(top, topdown=True):
    """
    os.path.walk like function for TDirectories.
    Return 4-tuple: (dirpath, dirnames, filenames, top)
        dirpath = 'file_name.root:/some/path' # may end in a '/'?
        dirnames = ['list', 'of' 'TDirectory', 'keys']
        filenames = ['list', 'of' 'object', 'keys']
        top = this level's TDirectory
    """
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    assert isinstance(top, rt.TDirectory)
    names = [k.GetName() for k in top.GetListOfKeys()]
    dirpath = top.GetPath()
    dirnames = []
    filenames = []
    ## filter names for directories
    for k in names:
        d = top.Get(k)
        if isinstance(d, rt.TDirectory):
            dirnames.append(k)
        else:
            filenames.append(k)
    ## sort
    dirnames.sort()
    filenames.sort()
    ## yield
    if topdown:
        yield dirpath, dirnames, filenames, top
    for dn in dirnames:
        d = top.Get(dn)
        for x in walk(d, topdown):
            yield x
    if not topdown:
        yield dirpath, dirnames, filenames, top

    
def convertTree2Dict(tree):
    
    runLumiDict = {}

    if not (hasattr(tree,'lumi') and hasattr(tree,'run')):
        print "tree does not contain run and lumi branches, returning empty json"
        return runLumiDict
    # loop over tree to get run, lumi "flat" dictionary
    tree.Draw('>>elist','','entrylist')        
    elist = rt.gDirectory.Get('elist')    
    entry = -1;
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)
        if entry%10000==0:
            print "processing entry %i"%entry
        if '%s'%(tree.run) in runLumiDict.keys():
            currentLumi = runLumiDict['%s'%(tree.run)]
            if int(tree.lumi) in currentLumi:
                pass
            else:                
                currentLumi.append(int(tree.lumi))
                runLumiDict.update({'%s'%(tree.run):currentLumi})
        else:
            runLumiDict['%s'%(tree.run)] = [int(tree.lumi)]

    # fix run, lumi list by grouping consecutive lumis
    for run in runLumiDict.keys():
        lumiGroups = []
        for k, g in groupby(enumerate(runLumiDict[run]), lambda (i,x):i-x):
            consecutiveLumis = map(itemgetter(1), g)
            lumiGroups.append([consecutiveLumis[0],consecutiveLumis[-1]])
        runLumiDict.update({run:lumiGroups})
        
    return runLumiDict
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-o','--output',dest="output",type="string",default="test.json",
                  help="Name of the json file to write to")
    
    (options,args) = parser.parse_args()

    
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)

    trees = []
    # crawl root file to look for trees
    for dirpath, dirnames, filenames, tdirectory in walk(rootFile):
        for filename in filenames:
            obj = tdirectory.Get(filename)
            if isinstance(obj, rt.TTree):
                trees.append(obj)
                

    # use first tree found
    runLumiDict = convertTree2Dict(trees[0])
    output = open(options.output,'w')
    json.dump(runLumiDict,output,sort_keys=True)
    output.close()
    print '\njson:'
    os.system('cat %s'%options.output)
    print '\n'
            
    
