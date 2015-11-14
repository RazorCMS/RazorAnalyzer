from optparse import OptionParser
import ROOT as rt
import sys
from array import *

#local imports
import rootTools
from framework import Config
from DustinTuple2RooDataSet import initializeWorkspace
from DustinTuples2DataCard import convertTree2TH1
from RunCombine import exec_me

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store datasets")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('--lumi-in',dest="lumi_in", default=1.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal systematic templates")
    parser.add_option('--num-pdf-weights',dest="numPdfWeights",default=60,type="int",
                  help="Number of nuisance parameters to use for PDF uncertainties")
    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    boxList = box.split('_')
    lumi = options.lumi
    lumi_in = options.lumi_in
    f = args[0]
    print 'Input file is %s' % f

    if options.noSignalSys:
        shapes = []
    else:
        shapes = ['muoneff','eleeff','jes','muontrig','eletrig','btag','muonfastsim','elefastsim','btagfastsim','facscale','renscale','facrenscale','ees','mes','lumi']
        shapes.extend([str(n)+'pdf' for n in range(options.numPdfWeights)])

    for curBox in boxList:
        #create workspace
        w = rt.RooWorkspace("w"+curBox)
        variables = initializeWorkspace(w,cfg,curBox)    
        btagMin =  w.var('nBtag').getMin()
        btagMax =  w.var('nBtag').getMax()
        z = array('d', cfg.getBinning(curBox)[2]) # nBtag binning
        
        #list of histograms
        ds = []
            
        #make MC histograms
        model = ''
        if f.lower().endswith('.root'):
            rootFile = rt.TFile.Open(f) #open file
            tree = rootFile.Get('RazorInclusive') #get tree

            #get histograms for sum of pdf and scale weights
            if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes or '0pdf' in shapes:
                nevents = rootFile.Get('NEvents')
                assert nevents
            else:
                nevents = None

            if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes:
                sumScaleWeights = rootFile.Get('SumScaleWeights')
                assert sumScaleWeights
            else:
                sumScaleWeights = None

            if '0pdf' in shapes:
                sumPdfWeights = rootFile.Get('SumPdfWeights')
                assert sumPdfWeights
            else:
                sumPdfWeights = None

            # get mass point information
            modelString = f.split('/')[-1].split('.root')[0].split('_')[0]
            model = modelString.split('-')[-1]
            massPoint = '_'.join(f.split('/')[-1].split('.root')[0].split('_')[1:3])
                               
            thyXsec = -1
            thyXsecErr = -1
            mGluino = -1
            mStop = -1
            mLSP = massPoint.split("_")[1]
            if "T1" in model:
                mGluino = massPoint.split("_")[0]
            if "T2" in model:
                mStop = massPoint.split("_")[0]
    
            if mGluino!=-1:
                for line in open('data/gluino13TeV.txt','r'):
                    line = line.replace('\n','')
                    if str(int(mGluino))==line.split(',')[0]:
                        thyXsec = float(line.split(',')[1]) #pb
                        thyXsecErr = 0.01*float(line.split(',')[2])
            if mStop!=-1:
                for line in open('data/stop13TeV.txt','r'):
                    line = line.replace('\n','')
                    if str(int(mStop))==line.split(',')[0]:
                        thyXsec = float(line.split(',')[1]) #pb
                        thyXsecErr = 0.01*float(line.split(',')[2]) 

            if isinstance( rootFile.Get('NEvents'), rt.TH1 ):
                nEvents = rootFile.Get('NEvents').Integral()
                globalScaleFactor = thyXsec*lumi/lumi_in/nEvents # FastSim samples
            else:
                globalScaleFactor = lumi/lumi_in # FullSim samples
                
            #get gluino and LSP masses
            tree.GetEntry(0)

            #add histogram to output file
            print("Building histogram for "+model)
            ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=globalScaleFactor, treeName=curBox+"_"+model))
            for shape in shapes:
                for updown in ["Up", "Down"]:
                    print("Building histogram for "+model+"_"+shape+updown)
                    ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=globalScaleFactor, treeName=curBox+"_"+model+"_"+shape+updown, sysErrOpt=shape+updown, sumScaleWeights=sumScaleWeights, sumPdfWeights=sumPdfWeights, nevents=nevents))
            rootFile.Close()
        else:
            print "Error: expected ROOT file!"
            sys.exit()

        #output file name
        if btagMax>btagMin+1:
            if "T1" in modelString:
                outFileName = '%s_%i_%i_lumi-%.3f_%i-%ibtag_%s.root'%(modelString,int(mGluino),int(mLSP),lumi/1000.,btagMin,btagMax-1,curBox)
            else:
                outFileName = '%s_%i_%i_lumi-%.3f_%i-%ibtag_%s.root'%(modelString,int(mStop),int(mLSP),lumi/1000.,btagMin,btagMax-1,curBox)
        else:
            if "T1" in modelString:
                outFileName = '%s_%i_%i_lumi-%.3f_%ibtag_%s.root'%(modelString,int(mGluino),int(mLSP),lumi/1000.,btagMin,curBox)
            else:
                outFileName = '%s_%i_%i_lumi-%.3f_%ibtag_%s.root'%(modelString,int(mStop),int(mLSP),lumi/1000.,btagMin,curBox)

        #output file
        print "Output File: %s"%(options.outDir+"/"+outFileName)
        outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
        outFile.cd()

        for hist in ds:
            print("Writing histogram: "+hist.GetName())
            hist.Write()
       
        outFile.Close()
