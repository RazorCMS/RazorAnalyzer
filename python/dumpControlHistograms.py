### INPUT: ROOT file produced by RazorAnalyzer python macro code

import sys
import argparse
import ROOT as rt

from macro import macro

parser = argparse.ArgumentParser()
parser.add_argument("fname",help="input ROOT file name")
parser.add_argument("varname", help="variable to print")
parser.add_argument("--all", dest="print_all",
        help="Print individual MC histograms")
args = parser.parse_args()
hists = macro.importHists(args.fname)

tot = hists["Data"][args.varname].Clone("tot")
tot.Reset()
for proc in hists:
    if args.print_all or proc == "Data":
        print proc
        hists[proc][args.varname].Print("all")
    else:
        tot.Add(hists[proc][args.varname])
print "Total MC"
tot.Print("all")
