from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config

k_T = 689.1/424.5
k_Z = 3.*2008.4/5482.
k_W = 3.*20508.9/50100.0
k_QCD = {'MultiJet': 1.2,
         'LooseLeptonMultiJet': 1.0,
         'DiJet': 1.2,
         'EleMultiJet':1.0,
         'MuMultiJet':1.0,
         'EleJet':1.0,
         'MuJet':1.0,
         'MuEle':1.0,
         'EleEle':1.0,
         'MuMu':1.0}
    
boxes = {'MuEle':0,
         'MuMu':1,
         'EleEle':2,
         'MuMultiJet':3,
         'MuJet':4,
         'EleMultiJet':5,
         'EleJet':6,
         'LooseLeptonMultiJet':7,
         'MultiJet':8,
         'DiJet':9}

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

def convertTree2Dataset(tree, cfg, box, workspace, useWeight, f, lumi, lumi_in, treeName='RMRTree'):
    """This defines the format of the RooDataSet"""
    
    if f.find('TTJets')!=-1:
        k = k_T*k_QCD[box]
    elif f.find('DYJets')!=-1 or f.find('ZJets')!=-1:
        k = k_Z*k_QCD[box]
    elif f.find('WJets')!=-1:
        k = k_W*k_QCD[box]
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
    

    btagCutoff = 3
    if box in ["MuEle", "MuMu", "EleEle"]:
        btagCutoff = 1

    tree.Draw('>>elist',
              'MR > %f && MR < %f && Rsq > %f && Rsq < %f && nBTaggedJets >= %f && box == %i && abs(dPhiRazor) < %f' % (mRmin,mRmax,rsqMin,rsqMax,btagMin,boxes[box],dPhiCut),
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
        if tree.nBTaggedJets >= btagCutoff:
            a.setRealValue('nBtag',btagCutoff)
        else:
            a.setRealValue('nBtag',tree.nBTaggedJets)
        if useWeight:
            try:
                a.setRealValue('W',tree.weight*lumi*k/lumi_in)
            except AttributeError:
                a.setRealValue('W',1.0)
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
    parser.add_option('-l','--lumi',dest="lumi", default=4000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--lumi-in',dest="lumi_in", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    lumi = options.lumi
    lumi_in = options.lumi_in
    useWeight = options.useWeight
    
    print 'Input files are %s' % ', '.join(args)
    
    w = rt.RooWorkspace("w"+box)

    initializeWorkspace(w,cfg)
    
    
    w.factory('W[1.,0.,+INF]')
    w.set('variables').add(w.var('W'))

    
    ds = []

    i = 0    
    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            tree = rootFile.Get('RazorInclusive')
            if f.lower().find('sms')==-1:
                i+=1
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
            
    if len(inFiles)==1:
        outFile = inFiles[0].split('/')[-1].replace('.root','_lumi-%.1f_%s.root'%(lumi/1000.,box))
    else:
        outFile = 'RazorAnalysis_SMCocktail_weighted_lumi-%.1f_%s.root'%(lumi/1000.,box)
        

    outFile = rt.TFile.Open(options.outDir+"/"+outFile,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
    
