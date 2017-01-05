import argparse
import ROOT as rt

from macro.razorFits import FitInstance

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('box', help='Analysis box')
    parser.add_argument('--tag', help='Analysis tag', 
            default='Razor2016_MoriondRereco')
    parser.add_argument('--mc', help='fit MC',
            action='store_true')
    parser.add_argument('--full', help='full fit',
            action='store_true')
    parser.add_argument('--load', help='load dataset from file',
            action='store_true')
    parser.add_argument('--load-fit', dest='loadFit', action='store_true',
            help='Create dataset but do not fit')
    parser.add_argument('--no-plot', dest='noPlot', action='store_true',
            help='Do not make plots')
    parser.add_argument('--unblind', action='store_true',
            help='Do not blind signal sensitive region')
    parser.add_argument('--run-toys', action='store_true',
            help='Run toy generation', dest='runToys')
    parser.add_argument('--plot-toys', action='store_true',
            help='Plot uncertainties from toys', dest='plotToys')
    args = parser.parse_args()

    weights = {}
    fitter = FitInstance(args.box, tag=args.tag, isData=not args.mc, 
            weights=weights, full=args.full)
    fitter.doFitSequence(load=args.load, doFit=(not args.loadFit),
            plot=(not args.noPlot), unblind=args.unblind, runToys=args.runToys,
            plotToys=args.plotToys)
