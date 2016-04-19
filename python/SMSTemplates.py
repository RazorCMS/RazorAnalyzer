from optparse import OptionParser
import ROOT as rt
import sys
from array import *

#local imports
import rootTools
from framework import Config
from DustinTuple2RooDataSet import initializeWorkspace
from DustinTuples2DataCard import convertTree2TH1, uncorrelate
from RunCombine import exec_me
from macro.razorAnalysis import xbinsSignal, colsSignal

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
    parser.add_option('--xBR',dest="xBR", default=-1,type="float",
                  help="x = BR(~g -> b b ~chi0)")
    parser.add_option('--yBR',dest="yBR", default=-1,type="float",
                  help="y = BR(~g -> t t ~chi0)")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('--no-signal-sys',dest="noSignalSys",default=False,action='store_true',
                  help="no signal systematic templates")
    parser.add_option('--merge-bins',dest="mergeBins", action="store_true",
                  help="merge some bins in Rsq")
    #pdf uncertainty options.  current prescription is just to take 10% uncorrelated error on each bin
    #parser.add_option('--num-pdf-weights',dest="numPdfWeights",default=0,type="int",
                  #help="Number of nuisance parameters to use for PDF uncertainties")
    #parser.add_option('--compute-pdf-envelope',dest="computePdfEnvelope",default=False,action='store_true',
                  #help="Use the SUS pdf reweighting prescription, summing weights in quadrature")
    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)

    box =  options.box
    boxList = box.split('_')
    lumi = options.lumi
    lumi_in = options.lumi_in
    f = args[0]
    print 'Input file is %s' % f

    unrollBins = None
    if options.mergeBins:
        btagBins = cfg.getBinning(box)[2][:-1]
        unrollBins = [(xbinsSignal[box][str(int(btags))+'B'], colsSignal[box][str(int(btags))+'B']) for btags in btagBins]

    if options.noSignalSys:
        shapes = []
    else:
        shapes = ['tightmuoneff','tighteleeff','vetomuoneff','vetoeleeff','jes','muontrig','eletrig','btag','tightmuonfastsim','tightelefastsim','vetomuonfastsim','vetoelefastsim','btagfastsim','facscale','renscale','facrenscale','ees','mes','pileup','isr','mcstat%s'%box.lower()]
        #shapes.append('pdf%s'%box.lower()) #this is for the flat 10% PDF uncertainty (uncorrelated across bins)
        #shapes.extend(['n'+str(n)+'pdf' for n in range(options.numPdfWeights)])

    for curBox in boxList:
        #create workspace
        w = rt.RooWorkspace("w"+curBox)
        variables = initializeWorkspace(w,cfg,curBox)    
        btagMin =  w.var('nBtag').getMin()
        btagMax =  w.var('nBtag').getMax()
        z = array('d', cfg.getBinning(curBox)[2]) # nBtag binning
        
        #list of histograms
        ds = []
        #dictionary of histograms (same content - useful for uncorrelate function)
        dsDict = {}
            
        #make MC histograms
        model = ''
        if f.lower().endswith('.root'):
            rootFile = rt.TFile.Open(f) #open file
            tree = rootFile.Get('RazorInclusive') #get tree

            #get histograms for sum of pdf and scale weights
            if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes or 'n0pdf' in shapes:
                nevents = rootFile.Get('NEvents')
                assert nevents
            else:
                nevents = None

            if 'facscale' in shapes or 'renscale' in shapes or 'facrenscale' in shapes:
                sumScaleWeights = rootFile.Get('SumScaleWeights')
                assert sumScaleWeights
            else:
                sumScaleWeights = None

            if 'n0pdf' in shapes:
                sumPdfWeights = rootFile.Get('SumPdfWeights')
                assert sumPdfWeights
            else:
                sumPdfWeights = None

            # get mass point information
            modelString = '_'.join(f.split('/')[-1].split('.root')[0].split('_')[:-2])
            if options.xBR>-1 and options.yBR>-1:
                modelString = modelString.replace('T1ttbb',('T1x%.2fy%.2f'%(options.xBR,options.yBR)).replace('.','p'))
            model = modelString.split('-')[-1]
            massPoint = '_'.join(f.split('/')[-1].split('.root')[0].split('_')[1:])
                               
            thyXsec = -1
            thyXsecErr = -1
            mGluino = -1
            mStop = -1
            mLSP = massPoint.split("_")[-1]
            if "T1" in model or "T5" in model:
                mGluino = massPoint.split("_")[-2]
            elif "T2" in model:
                mStop = massPoint.split("_")[-2]
    
            if mGluino!=-1:
                for line in open('data/gluino13TeV.txt','r'):
                    line = line.replace('\n','')
                    if str(int(mGluino))==line.split(',')[0]:
                        thyXsec = float(line.split(',')[1]) #pb
                        thyXsecErr = 0.01*float(line.split(',')[2])
                if model=="T5ttttDM175T2tt":                    
                    for line in open('data/stop13TeV.txt','r'):
                        line = line.replace('\n','')
                        if str(int(mLSP)+175)==line.split(',')[0]:
                            thyXsecStop = float(line.split(',')[1]) #pb
                    thyXsec+=thyXsecStop                    
            elif mStop!=-1:
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
            ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=globalScaleFactor, treeName=curBox+"_"+model, unrollBins=unrollBins, xBR=options.xBR, yBR=options.yBR))
            for shape in shapes:
                for updown in ["Up", "Down"]:
                    print("Building histogram for "+model+"_"+shape+updown)
                    ds.append(convertTree2TH1(tree, cfg, curBox, w, f, globalScaleFactor=globalScaleFactor, treeName=curBox+"_"+model+"_"+shape+updown, sysErrOpt=shape+updown, sumScaleWeights=sumScaleWeights, sumPdfWeights=sumPdfWeights, nevents=nevents, unrollBins=unrollBins, xBR=options.xBR, yBR=options.yBR))
            rootFile.Close()

            #make pdf envelope up/down (for SUS pdf uncertainty prescription)
            #if options.computePdfEnvelope:
            #    print "Building pdf envelope up/down histograms"
            #    #make summed histogram
            #    pdfEnvelopeUp = ds[0].Clone(ds[0].GetName()+'_pdfenvelopeUp')
            #    pdfEnvelopeDown = ds[0].Clone(ds[0].GetName()+'_pdfenvelopeDown')
            #    pdfEnvelopeUp.Reset()
            #    pdfEnvelopeDown.Reset()
            #    for p in range(options.numPdfWeights):
            #        print "Adding pdf variation",p,"to envelope"
            #        #get correct pdf histogram
            #        thisVariationUp = None
            #        thisVariationDown = None
            #        for h in ds:
            #            if 'n'+str(p)+'pdfUp' in h.GetName():
            #                thisVariationUp = h
            #            elif 'n'+str(p)+'pdfDown' in h.GetName():
            #                thisVariationDown = h
            #        if thisVariationUp is None or thisVariationDown is None:
            #            print "Error: did not find pdf variation histogram",p
            #            continue
            #        for bx in range(1,ds[0].GetNbinsX()+1):
            #            #add (up - nominal) in quadrature
            #            newUp = ((pdfEnvelopeUp.GetBinContent(bx))**2 + (thisVariationUp.GetBinContent(bx) - ds[0].GetBinContent(bx))**2)**(0.5)   
            #            pdfEnvelopeUp.SetBinContent(bx, newUp)
            #            #add (down - nominal) in quadrature
            #            newDown = -((pdfEnvelopeDown.GetBinContent(bx))**2 + (thisVariationDown.GetBinContent(bx) - ds[0].GetBinContent(bx))**2)**(0.5)   
            #            pdfEnvelopeDown.SetBinContent(bx, newDown)
            #    pdfEnvelopeUp.Add(ds[0])
            #    pdfEnvelopeDown.Add(ds[0])
            #    #zero any negative bins
            #    for bx in range(1,ds[0].GetNbinsX()+1):
            #        if pdfEnvelopeDown.GetBinContent(bx) < 0:
            #            pdfEnvelopeDown.SetBinContent(bx,0)
            #    #normalize
            #    pdfEnvelopeUp.Scale(ds[0].Integral()*1.0/pdfEnvelopeUp.Integral())
            #    pdfEnvelopeDown.Scale(ds[0].Integral()*1.0/pdfEnvelopeDown.Integral())
            #    #append
            #    ds.append(pdfEnvelopeUp)
            #    ds.append(pdfEnvelopeDown)
                
            #convert dataset list to dict
            for d in ds: dsDict[d.GetName()] = d
                
            #perform uncorrelation procedure (for MC stat and pdf uncertainty)
            if 'pdf%s'%box.lower() in shapes:
                uncorrelate(dsDict, 'pdf%s'%box.lower())
                # remove empty bins 
                for bx in range(1,ds[0].GetNbinsX()+1):
                    if ds[0].GetBinContent(bx) == 0 or ds[0].GetBinError(bx) == 0:
                        print 'Removing empty PDF uncertainty bin %i'%(bx)
                        del dsDict['%s_pdf%s%iUp'%(ds[0].GetName(),box.lower(),bx)]
                        del dsDict['%s_pdf%s%iDown'%(ds[0].GetName(),box.lower(),bx)]
            if 'mcstat%s'%box.lower() in shapes:
                uncorrelate(dsDict, 'mcstat%s'%box.lower())
                # remove unnecessary MC stat bins (relative uncertainty < 10%) see htt recommendation
                # https://indico.cern.ch/event/373752/session/6/contribution/14/attachments/744534/1021298/bbb-HCG.pdf
                for bx in range(1,ds[0].GetNbinsX()+1):
                    if ds[0].GetBinContent(bx) == 0 or ds[0].GetBinError(bx) == 0:
                        print 'Relative MC stat uncertainty bin %i = %.1f%% is less than 10%%'%(bx,0.)
                        print 'Removing histogram: %s_mcstat%s%iUp'%(ds[0].GetName(),box.lower(),bx)
                        print 'Removing histogram: %s_mcstat%s%iDown'%(ds[0].GetName(),box.lower(),bx)
                        del dsDict['%s_mcstat%s%iUp'%(ds[0].GetName(),box.lower(),bx)]
                        del dsDict['%s_mcstat%s%iDown'%(ds[0].GetName(),box.lower(),bx)]
                    elif ds[0].GetBinContent(bx) > 0 and ds[0].GetBinError(bx)/ds[0].GetBinContent(bx) < 0.1:
                        print 'Relative MC stat uncertainty bin %i = %.1f%% is less than 10%%'%(bx,100.*ds[0].GetBinError(bx)/ds[0].GetBinContent(bx))
                        print 'Removing histogram: %s_mcstat%s%iUp'%(ds[0].GetName(),box.lower(),bx)
                        print 'Removing histogram: %s_mcstat%s%iDown'%(ds[0].GetName(),box.lower(),bx)
                        del dsDict['%s_mcstat%s%iUp'%(ds[0].GetName(),box.lower(),bx)]
                        del dsDict['%s_mcstat%s%iDown'%(ds[0].GetName(),box.lower(),bx)]
        else:
            print "Error: expected ROOT file!"
            sys.exit()

        #output file name
        if btagMax>btagMin+1:
            if "T1" in modelString or "T5" in modelString:
                outFileName = '%s_%i_%i_lumi-%.3f_%i-%ibtag_%s.root'%(modelString,int(mGluino),int(mLSP),lumi/1000.,btagMin,btagMax-1,curBox)
            else:
                outFileName = '%s_%i_%i_lumi-%.3f_%i-%ibtag_%s.root'%(modelString,int(mStop),int(mLSP),lumi/1000.,btagMin,btagMax-1,curBox)
        else:
            if "T1" in modelString or "T5" in modelString:
                outFileName = '%s_%i_%i_lumi-%.3f_%ibtag_%s.root'%(modelString,int(mGluino),int(mLSP),lumi/1000.,btagMin,curBox)
            else:
                outFileName = '%s_%i_%i_lumi-%.3f_%ibtag_%s.root'%(modelString,int(mStop),int(mLSP),lumi/1000.,btagMin,curBox)

        #output file
        print "Output File: %s"%(options.outDir+"/"+outFileName)
        outFile = rt.TFile.Open(options.outDir+"/"+outFileName,'recreate')
        outFile.cd()

        
        for name in dsDict:            
            #if 'pdf' in name and options.computePdfEnvelope and not 'envelope' in name:
            #    continue
            print("Writing histogram: "+dsDict[name].GetName())
            dsDict[name].Write()

       
        outFile.Close()
