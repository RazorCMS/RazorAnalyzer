import argparse
import ROOT as rt

from macro.razorFits import FitInstance

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Configuration
    parser.add_argument('box', help='Analysis box')
    parser.add_argument('--tag', help='Analysis tag', 
            default='Razor2016_MoriondRereco')
    parser.add_argument('--config', help='Fit config to use', 
            default='config/run2_2017_01_07_Run2016G_SUSYUnblind_Sep23ReReco.config')
    # Fit types
    parser.add_argument('--mc', help='fit MC',
            action='store_true')
    parser.add_argument('--full', help='full fit',
            action='store_true')
    # Actions
    parser.add_argument('--load', help='load dataset from file',
            action='store_true')
    parser.add_argument('--load-fit', dest='loadFit', action='store_true',
            help='Use existing fit result')
    parser.add_argument('--no-plot', dest='noPlot', action='store_true',
            help='Do not make plots')
    parser.add_argument('--unblind', action='store_true',
            help='Do not blind signal sensitive region')
    parser.add_argument('--run-toys', action='store_true',
            help='Run toy generation', dest='runToys')
    parser.add_argument('--load-toys', action='store_true',
            help='Plot uncertainties from toys', dest='loadToys')
    parser.add_argument('--input-fit-file', dest='inputFitFile',
            help='File to load existing fit from')
    parser.add_argument('--no-fit', dest='noFit', action='store_true',
            help='Do not perform fit')
    args = parser.parse_args()

    weights = {}
    fitter = FitInstance(args.box, tag=args.tag, isData=not args.mc, 
            weights=weights, full=args.full, configFile=args.config)
    fitter.doFitSequence(load=(args.load or args.loadFit), 
            doFit=(not (args.loadFit or args.noFit)), plot=(not args.noPlot), 
            unblind=args.unblind, runToys=args.runToys, loadToys=args.loadToys, 
            inputFitFile=args.inputFitFile)
