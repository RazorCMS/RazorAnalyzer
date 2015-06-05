from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
import sys
from array import *

k_T = 689.1/424.5
k_Z = 3.*2008.4/5482.
k_W = 3.*20508.9/50100.0

k_QCD = {}
    
boxes = {'MuEle':[0],
         'MuMu':[1],
         'EleEle':[2],
         'MuSixJet':[3],
         'MuFourJet':[4],
         'MuMultiJet':[3,4],
         'MuJet':[5],
         'EleSixJet':[6],
         'EleFourJet':[7],
         'EleMultiJet':[6,7],
         'EleJet':[8],
         'LooseLeptonSixJet':[9],
         'LooseLeptonFourJet':[10],
         'LooseLeptonMultiJet':[9,10],
         'SixJet':[11],
         'FourJet':[12],
         'MultiJet':[11,12],
         'LooseLeptonDiJet':[13],
         'DiJet':[14]}

dPhiCut = 2.7

def initializeWorkspace(w,cfg):
    variables = cfg.getVariablesRange(box,"variables",w)
    parameters = cfg.getVariables(box, "parameters")
    paramNames = []
    for parameter in parameters:
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if paramName.find("Cut")==-1 and paramName.find("Ntot")==-1:
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
        else:
            if paramName.find("Ntot")==-1:
                w.var(paramName).setConstant(True)
            else:
                w.var(paramName).setConstant(False)


def getSumOfWeights(tree, cfg, box, workspace, useWeight, f, lumi, lumi_in):
    if f.find('SMS')!=-1:
        k = 1.
    elif f.find('TTJets')!=-1:
        k = k_T
    elif f.find('DYJets')!=-1 or f.find('ZJets')!=-1:
        k = k_Z
    elif f.find('WJets')!=-1:
        k = k_W
    else:
        k = 1.
        
    args = workspace.set("variables")
    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin()
    mRmax = args['MR'].getMax()

    #we cut away events outside our Rsq window
    rsqMin = args['Rsq'].getMin()
    rsqMax = args['Rsq'].getMax()

    btagMin =  args['nBtag'].getMin()
    btagMax =  args['nBtag'].getMax()
    
    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    btagCutoff = 3
    if box in ["MuEle", "MuMu", "EleEle"]:
        btagCutoff = 1
        
    boxCut = '(' + ' || '.join(['box==%i'%boxNum for boxNum in boxes[box]]) + ')'

    label = f.replace('.root','').split('/')[-1]
    htemp = rt.TH1D('htemp_%s'%label,'htemp_%s'%label,len(z)-1,z)

    if useWeight:
        tree.Project(htemp.GetName(),
                    'min(nBTaggedJets,%i)'%btagCutoff,
                    '(%f/%f) * %f * weight * (MR > %f && MR < %f && Rsq > %f && Rsq < %f && min(nBTaggedJets,%i) >= %i && min(nBTaggedJets,%i) < %f && %s && abs(dPhiRazor) < %f)' % (lumi,lumi_in,k,mRmin,mRmax,rsqMin,rsqMax,btagCutoff,btagMin,btagCutoff,btagMax,boxCut,dPhiCut))
    else:
        tree.Project(htemp.GetName(),
                    'MR',
                    '(MR > %f && MR < %f && Rsq > %f && Rsq < %f && min(nBTaggedJets,%i) >= %i && min(nBTaggedJets,%i) < %f && %s && abs(dPhiRazor) < %f)' % (mRmin,mRmax,rsqMin,rsqMax,btagCutoff,btagMin,btagCutoff,btagMax,boxCut,dPhiCut))
        
    return [htemp.GetBinContent(i) for i in range(1,len(z))]
        
    
def convertTree2Dataset(tree, cfg, box, workspace, useWeight, f, lumi, lumi_in, treeName='RMRTree'):
    """This defines the format of the RooDataSet"""
    if f.find('SMS')!=-1:
        k = [1. for k_btag in k_QCD[box]]
    elif f.find('TTJets')!=-1:
        k = [k_T*k_btag for k_btag in k_QCD[box]]
    elif f.find('DYJets')!=-1 or f.find('ZJets')!=-1:
        k = [k_Z*k_btag for k_btag in k_QCD[box]]
    elif f.find('WJets')!=-1:
        k = [k_W*k_btag for k_btag in k_QCD[box]]
    else:
        k = k_QCD[box]

    args = workspace.set("variables")
    data = rt.RooDataSet(treeName,'Selected R and MR',args)
    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin()
    mRmax = args['MR'].getMax()

    #we cut away events outside our Rsq window
    rsqMin = args['Rsq'].getMin()
    rsqMax = args['Rsq'].getMax()

    btagMin =  args['nBtag'].getMin()
    btagMax =  args['nBtag'].getMax()

    z = array('d', cfg.getBinning(box)[2]) # nBtag binning
    
    label = f.replace('.root','').split('/')[-1]
    htemp = rt.TH1D('htemp2_%s'%label,'htemp2_%s'%label,len(z)-1,z)

    btagCutoff = 3
    if box in ["MuEle", "MuMu", "EleEle"]:
        btagCutoff = 1
        
    boxCut = '(' + ' || '.join(['box==%i'%boxNum for boxNum in boxes[box]]) + ')'

    tree.Draw('>>elist',
              'MR > %f && MR < %f && Rsq > %f && Rsq < %f && min(nBTaggedJets,%i) >= %i && min(nBTaggedJets,%i) < %i && %s && abs(dPhiRazor) < %f' % (mRmin,mRmax,rsqMin,rsqMax,btagCutoff,btagMin,btagCutoff,btagMax,boxCut,dPhiCut),
              'entrylist')
        
    elist = rt.gDirectory.Get('elist')
    
    entry = -1;
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)

        #set the RooArgSet and save
        a = rt.RooArgSet(args)
        
        a.setRealValue('MR',tree.MR)
        a.setRealValue('Rsq',tree.Rsq)
        a.setRealValue('nBtag',min(tree.nBTaggedJets,btagCutoff))
        
        if useWeight:
            btag_bin = htemp.FindBin(min(tree.nBTaggedJets,btagCutoff)) - 1
            a.setRealValue('W',tree.weight*lumi*k[btag_bin]/lumi_in)
        else:
            a.setRealValue('W',1.0)
        data.add(a)
        
    numEntries = data.numEntries()
    
    wdata = rt.RooDataSet(data.GetName(),data.GetTitle(),data,data.get(),"MR>=0.","W")
    print "Filename: %s"%f
    print "Number of Entries in Box %s = %d"%(box,data.numEntries())
    print "Sum of Weights in Box %s = %.1f"%(box,wdata.sumEntries())

    return wdata

    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-w','--weight',dest="useWeight",default=False,action='store_true',
                  help="use weight")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--lumi-in',dest="lumi_in", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-q','--remove-qcd',dest="removeQCD",default=False,action='store_true',
                  help="remove QCD, while augmenting remaining MC backgrounds")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    lumi = options.lumi
    lumi_in = options.lumi_in
    useWeight = options.useWeight
    removeQCD = options.removeQCD
    
    print 'Input files are %s' % ', '.join(args)
    
    w = rt.RooWorkspace("w"+box)

    initializeWorkspace(w,cfg)
    
    
    w.factory('W[1.,0.,+INF]')
    w.set('variables').add(w.var('W'))
    
    ds = []

    
    btagMin =  w.var('nBtag').getMin()
    btagMax =  w.var('nBtag').getMax()

    if removeQCD:
        # first get sum of weights for each background per b-tag bin ( sumW[label] )
        sumW = {}
        sumWQCD = 0.
        for f in args:
            if f.lower().endswith('.root'):
                rootFile = rt.TFile(f)
                tree = rootFile.Get('RazorInclusive')
                if f.lower().find('sms')==-1:
                    
                    label = f.replace('.root','').split('/')[-1]
                    sumW[label] = getSumOfWeights(tree, cfg, box, w, useWeight, f, lumi, lumi_in)
                    if label.find('QCD')!=-1: sumWQCD = sumW[label]
        # get total sum of weights
        sumWTotal = [sum(allW) for allW in zip( * sumW.values() )]

        # get scale factor to scale other backgrounds by
        k_QCD[box] = [total/(total - qcd) for total, qcd in zip(sumWTotal,sumWQCD)]
         
        print "sum of weights total     =", sumWTotal
        print "sum of weights QCD       =", sumWQCD
        print "scale factor k_QCD[ %s ] ="%box, k_QCD[box]

    for i, f in enumerate(args):
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            tree = rootFile.Get('RazorInclusive')
            if f.lower().find('sms')==-1:
                if removeQCD and f.find('QCD')!=-1:
                    continue # do not add QCD
                else:
                    ds.append(convertTree2Dataset(tree, cfg, box, w, useWeight, f, lumi, lumi_in,  'RMRTree_%i'%i))
                
            else:
                model = f.split('-')[1].split('_')[0]
                massPoint = '_'.join(f.split('_')[2:4])
                ds.append(convertTree2Dataset(tree, cfg, box, w, useWeight, f ,lumi, lumi_in, 'signal'))
                
    wdata = ds[0].Clone('RMRTree')
    for ids in range(1,len(ds)):
        wdata.append(ds[ids])
    
    rootTools.Utils.importToWS(w,wdata)
    
    inFiles = [f for f in args if f.lower().endswith('.root')]
    
    args = w.set("variables")
    
    #we cut away events outside our MR window
    mRmin = args['MR'].getMin()
    mRmax = args['MR'].getMax()

    #we cut away events outside our Rsq window
    rsqMin = args['Rsq'].getMin()
    rsqMax = args['Rsq'].getMax()

    btagMin =  args['nBtag'].getMin()
    btagMax =  args['nBtag'].getMax()
    
    if len(inFiles)==1:
        if btagMax>btagMin+1:
            outFile = inFiles[0].split('/')[-1].replace('.root','_lumi-%.1f_%i-%ibtag_%s.root'%(lumi/1000.,btagMin,btagMax-1,box))
        else:
            outFile = inFiles[0].split('/')[-1].replace('.root','_lumi-%.1f_%ibtag_%s.root'%(lumi/1000.,btagMin,box))
    else:
        if btagMax>btagMin+1:
            outFile = 'RazorAnalysis_SMCocktail_weighted_lumi-%.1f_%i-%ibtag_%s.root'%(lumi/1000.,btagMin,btagMax-1,box)
        else:
            outFile = 'RazorAnalysis_SMCocktail_weighted_lumi-%.1f_%ibtag_%s.root'%(lumi/1000.,btagMin,box)
        

    print "Output file is: %s" % (options.outDir+"/"+outFile)
    outFile = rt.TFile.Open(options.outDir+"/"+outFile,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
    
