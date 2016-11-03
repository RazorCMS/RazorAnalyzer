import sys
import argparse
import ROOT as rt

from macro import macro

parser = argparse.ArgumentParser()
parser.add_argument("fname",help="input ROOT file name")
parser.add_argument("varname", help="variable to print")
args = parser.parse_args()
hists = macro.importHists(args.fname)

for proc in hists:
    print "Process",proc
    hists[proc][args.varname].Print("all")
