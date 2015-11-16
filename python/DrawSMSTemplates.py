import sys,os
import argparse
import copy
import ROOT as rt

#local imports
from framework import Config
import macro.plotting as plotting

if __name__ == "__main__":
    rt.gROOT.SetBatch()

    #parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('--num-pdf-weights',dest="numPdfWeights",default=60,type=int, help="Number of nuisance parameters to use for PDF uncertainties")
    parser.add_argument("filename", help="Path to input file")
    args = parser.parse_args() 

    #get input file
    infile = rt.TFile(args.filename)
    assert infile

    #specify shape uncertainties
    shapes = ['muoneff','eleeff','jes','muontrig','eletrig','btag','muonfastsim','elefastsim','btagfastsim','facscale','renscale','facrenscale','ees','mes']
    shapes.extend([str(n)+'pdf' for n in range(args.numPdfWeights)])
    titles = {'muoneff':'Muon efficiency', 'eleeff':'Electron efficiency', 'jes':'Jet energy scale',
            'muontrig':'Muon trigger eff.', 'eletrig':'Electron trigger eff.',
            'btag':'B-tagging efficiency', 'muonfastsim':'Fastsim muon eff.',
            'elefastsim':'Fastsim electron eff.', 'btagfastsim':'Fastsim b-tagging eff.',
            'facscale':'Factorization scale', 'renscale':'Renormalization scale',
            'facrenscale':'Diag. scale variation', 'ees':'Electron energy scale',
            'mes':'Muon energy scale'}
    for n in range(args.numPdfWeights):
        titles[str(n)+'pdf'] = 'Pdf variation '+str(n)

    #parse model and box name from file
    split = args.filename.replace('.root','').split('_')
    model = split[0].replace('SMS-','')
    box = split[-1]

    #make output directory
    dirName = "SMSPlots"
    os.system('mkdir -p '+dirName)

    #load histograms and plot
    hists = {}
    hists['nominal'] = infile.Get('_'.join([box,model]))
    print "Loading histogram",'_'.join([box,model])
    c = rt.TCanvas('c','c',800,600)
    assert hists['nominal']
    for shape in shapes:
        for updown in ['Up','Down']:
            hists[shape+updown] = infile.Get('_'.join([box,model,shape+updown]))
            print "Loading histogram",'_'.join([box,model,shape+updown])
            assert hists[shape+updown]
        toplot = [hists['nominal'],hists[shape+'Up'],hists[shape+'Down']]
        colors = [rt.kBlack, rt.kGreen, rt.kBlue]
        title = titles[shape]
        d = {0:toplot[0],1:toplot[1],2:toplot[2]}
        t = {0:'Nominal',1:title+' Up',2:title+' Down'}
        leg = plotting.makeLegend(d, t, [0,1,2], 0.4, 0.0, 0.6, 0.3)
        plotting.plot_several(c, toplot, leg, colors, 'Bin', 'A.U.', printstr=box+model+shape, lumistr='', commentstr=model+' '+box, printdir=dirName)
