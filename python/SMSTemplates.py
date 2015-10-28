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
    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    boxList = box.split('_')
    lumi = options.lumi
    lumi_in = options.lumi_in
    f = args[0]
    print 'Input file is %s' % f

    ##################
    #get pileup weight hist (remove this later)
    pileupWeightFileName = "root://eoscms:///eos/cms/store/group/phys_susy/razor/Run2Analysis/ScaleFactors/PileupWeights/NVtxReweight_ZToMuMu_2015D_1264ipb.root"
    pileupHistName = "NVtxReweight"
    pileupWeightFile = rt.TFile.Open(pileupWeightFileName)
    pileupWeightHist = pileupWeightFile.Get(pileupHistName)
    ##################

    #list of shape systematics to apply.
    shapes = ["muoneff", "eleeff", "jes"]

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
        modelString = ""
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f) #open file
            tree = rootFile.Get('RazorInclusive') #get tree
            
            #get gluino and LSP masses
            tree.GetEntry(0)
            modelString = (f.split('/')[-1]).replace('.root','') #assumes filename like SMS-T1bbbb_X_Y.root

            #add histogram to output file
            print("Building histogram for "+modelString)
            ds.append(convertTree2TH1(tree, cfg, curBox, w, f, lumi, lumi_in, modelString, pileupWeightHist=pileupWeightHist, hadronicTriggerWeight=0.935))
            for shape in shapes:
                for updown in ["Up", "Down"]:
                    print("Building histogram for "+modelString+"_"+shape+updown)
                    ds.append(convertTree2TH1(tree, cfg, curBox, w, f, lumi, lumi_in, modelString+"_"+shape+updown, sysErrOpt=shape+updown, pileupWeightHist=pileupWeightHist, hadronicTriggerWeight=0.935))
            rootFile.Close()
        else:
            print "Error: expected ROOT file!"
            sys.exit()

        #output file name
        if btagMax>btagMin+1:
            outFileName = 'RazorInclusive_%s_lumi-%.1f_%i-%ibtag_%s.root'%(modelString,lumi/1000.,btagMin,btagMax-1,curBox)
        else:
            outFileName = 'RazorInclusive_%s_lumi-%.1f_%ibtag_%s.root'%(modelString,lumi/1000.,btagMin,curBox)

        #output file
        print "Output File: %s"%(options.outDir+"/"+outFileName)
        outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
        outFile.cd()

        for hist in ds:
            print("Writing histogram: "+hist.GetName())
            hist.Write()
       
        outFile.Close()
