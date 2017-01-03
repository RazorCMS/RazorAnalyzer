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
    args = parser.parse_args()

    weights = { 
            #"ZInv":1.75 
            }
    fitter = FitInstance(args.box, tag=args.tag, isData=not args.mc, 
            weights=weights)
    fitter.initDataset()
    fitter.writeWorkspace()
