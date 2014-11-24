import sys, os
import string
import numpy as np
from scipy import stats
import itertools
import ROOT as rt
import rootTools
import math

def Gamma(a, x):
    return rt.TMath.Gamma(a) * rt.Math.inc_gamma_c(a,x)

def Gfun(x, y, X0, Y0, B, N):
    return Gamma(N,B*N*rt.TMath.Power((x-X0)*(y-Y0),1/N))

def RazorIntegral(B, N, X0, Y0, XMin, XMax, YMin, YMax):
    return (N/rt.TMath.Power(B*N,N))*(Gfun(XMin,YMin,X0,Y0,B,N)-Gfun(XMin,YMax,X0,Y0,B,N)-Gfun(XMax,YMin,X0,Y0,B,N)+Gfun(XMax,YMax,X0,Y0,B,N))

def convertDataset2Unweighted(wdata, box):
    """Get the cocktail dataset from the file"""
    
    row = wdata.get()

    MR = row['MR']
    Rsq = row['Rsq']
    nBTaggedJets = row['nBTaggedJets']
    
    varSet = rt.RooArgSet(MR,Rsq,nBTaggedJets)
    varList = rt.RooArgList(MR,Rsq,nBTaggedJets)
    varList2D = rt.RooArgList(MR,Rsq)
    uwdata = rt.RooDataSet('RMRTree','Unweighted Cocktail',varSet)
        
    mRmin = row['MR'].getMin()
    mRmax = row['MR'].getMax()
    rsqMin = row['Rsq'].getMin()
    rsqMax = row['Rsq'].getMax()
    nbtagMin = row['nBTaggedJets'].getMin()
    nbtagMax = row['nBTaggedJets'].getMax()
    
    myTH3 = rt.TH3D("h%s"%box, "h%s"%box, 100, 400, 4000, 70, 0.2, 1.5, 3, 1, 4)
    myTH2 = rt.TH2D("h", "h", 100, 400, 4000, 70, 0.2, 1.5)
    myTH2.Sumw2()

    # fills automatically with weight
    wdata.fillHistogram(myTH3, varList,"MR>0")
    wdata.fillHistogram(myTH2, varList2D,"MR>0")
    
    c = rt.TCanvas()
    rt.gStyle.SetOptStat(1001000011)
    myTH2.SetTitle("Weighted %s"%box)
    sumW2 = 0
    for i in range(0,wdata.numEntries()):
       wdata.get(i)
       sumW2+=(wdata.weight())*(wdata.weight())

    print "sum (weights)^2 = %.1f" %sumW2
    print "(sum weights)^2 = %.1f" %((wdata.sumEntries())*(wdata.sumEntries()))
    effEntries = (((wdata.sumEntries())*(wdata.sumEntries()))/sumW2)
    print "effective entries = %.1f"%effEntries
    myTH2.Draw("colz")
    c.SetLogy()
    c.SetLogx()
    c.Print("Cocktail_%s_DatasetWeighted.pdf"%box)

    print wdata.weight()
    Nev = myTH3.Integral()
    Nent = myTH3.GetEntries()
    print "weighted events %.1f"% Nev
    print "entries  %d"% Nent
    Npois = rt.RooRandom.randomGenerator().Poisson(Nev)
    for i in range(0,Npois):
       myMR = rt.Double()
       myRsq = rt.Double()
       mynBTaggedJets = rt.Double()
       myTH3.GetRandom3(myMR,myRsq,mynBTaggedJets)
       mynBTaggedJets = int(mynBTaggedJets)
       varSet.setRealValue('MR',myMR)
       varSet.setRealValue('Rsq',myRsq)
       varSet.setRealValue('nBTaggedJets',mynBTaggedJets)
       uwdata.add(varSet)

    myTH2Toy = rt.TH2D("h", "h", 100, 400, 4000, 70, 0.25, 1.5)
    uwdata.fillHistogram(myTH2Toy, varList2D,"MR>0")
    myTH2Toy.SetTitle("Unweighted %s"%box)
    myTH2Toy.Draw("colz")
    c.Print("Cocktail_%s_ToyUnweighted.pdf"%box)
   
    return uwdata

##initialization

#check that input file was specified
if len(sys.argv) < 2: 
    print 'Usage: RazorInclusiveFits.py filename1.root ...'
    exit()

#switch ROOT to batch mode
rt.gROOT.SetBatch()
#turn on fit stats
rt.gStyle.SetOptFit(0111)

#get the path of this script and load the razor library
pathname = os.path.dirname(sys.argv[0]) 
fullpath = os.path.abspath(pathname) #get the full path of this script
outpath = fullpath+"/output"
libpath = fullpath+"/lib"
rt.gSystem.Load(libpath+"/libRazorRun2.so")

#suppress info messages from RooFit
rt.RooMsgService.instance().setStreamStatus(1,False)

#load the TTrees from the input files
boxName = 'MultiJet'
print 'Will perform a fit in the '+boxName+' box and compute likelihoods for signal and background'
filenames = [sys.argv[i] for i in range(1, len(sys.argv))]
files = [rt.TFile(filename) for filename in filenames]
backgroundTrees = [file.Get(boxName) for index,file in enumerate(files) if 'SMS' not in filenames[index]]
signalTrees = [file.Get(boxName) for index,file in enumerate(files) if 'SMS' in filenames[index]]
if len(signalTrees) != 1: 
    print("Error: please include exactly one signal SMS file!")
    exit()
signalTree = signalTrees[0]

##set important parameters
intLuminosity = 5000 #in inverse picobarns
bTagBox = 3 #how many b-tags to require
sigma = 0.01 #in picobarns
nSig = 105964 #number of processed events in signal sample
signalWeight = sigma / nSig * intLuminosity
pOfSigma = 1.0 #prior probability for this value of sigma

##create RooDataSets with the MR and Rsq values of events in each tree

#set bins for MR and Rsq
MRMin = 400
RsqMin = 0.2
MRMax = 4000
RsqMax = 1.5
MRBinning = rt.RooBinning(MRMin, MRMax)
RsqBinning = rt.RooBinning(RsqMin, RsqMax)
MRBinLowEdges = [400, 450, 550, 700, 900, 1200, 1600, 2500]
MRBinHighEdges = [MRBinLowEdges[i] for i in range(1, len(MRBinLowEdges))]
MRBinHighEdges.append(MRMax)
RsqBinLowEdges = [0.2, 0.25, 0.30,0.41,0.52,0.64,0.80]
RsqBinHighEdges = [RsqBinLowEdges[i] for i in range(1, len(RsqBinLowEdges))]
RsqBinHighEdges.append(RsqMax)
for i in range(1, len(MRBinLowEdges)): MRBinning.addBoundary(MRBinLowEdges[i])
for i in range(1, len(RsqBinLowEdges)): RsqBinning.addBoundary(RsqBinLowEdges[i])
MR = rt.RooRealVar("MR", "MR", 400, MRMax) 
Rsq = rt.RooRealVar("Rsq", "Rsq", 0.2, RsqMax)
nBTaggedJets = rt.RooRealVar("nBTaggedJets", "nBTaggedJets", 0, 10)
MR.setBinning(MRBinning)
Rsq.setBinning(RsqBinning)

#define the fit region
MRSideband = 100
RsqSideband = 0.1

#create sideband and extrapolation regions
fitMRMax = 700
fitRsqMax = 0.6 
#sidebands (2 rectangles)
MR.setRange("lowMR", MRMin, MRMin+MRSideband)
Rsq.setRange("lowMR", RsqMin, fitRsqMax)
MR.setRange("lowRsq", MRMin, fitMRMax)
Rsq.setRange("lowRsq", RsqMin, RsqMin+RsqSideband)
#extrapolation region (3 rectangles)
MR.setRange("extr1", fitMRMax, MRMax) #region 1: high MR, low Rsq
Rsq.setRange("extr1", RsqMin, RsqMin+RsqSideband)
MR.setRange("extr2", MRMin, MRMin+MRSideband) #region 2: low MR, high Rsq
Rsq.setRange("extr2", fitRsqMax, RsqMax)
MR.setRange("extr3", MRMin+MRSideband, MRMax) #region 3: high MR, high Rsq
Rsq.setRange("extr3", RsqMin+RsqSideband, RsqMax)
#full region
MR.setRange("full", MRMin, MRMax)
Rsq.setRange("full", RsqMin, RsqMax)

binCorners = [[[MRBinMin,MRBinMax], [RsqBinMin, RsqBinMax]] for [(MRBinMin,MRBinMax), (RsqBinMin,RsqBinMax)] in itertools.product(zip(MRBinLowEdges,MRBinHighEdges), zip(RsqBinLowEdges,RsqBinHighEdges)) if MRBinMin+1 >= fitMRMax or RsqBinMin+0.01 >= fitRsqMax or (MRBinMin+1 >= MRMin+MRSideband and RsqBinMin+0.01 >= RsqMin+RsqSideband)] #list of all bins in the extrapolation region

#create RooDataSets for signal and background
weight = rt.RooRealVar("weight", "weight", 0.0, 100000.0)
args2D = rt.RooArgSet(MR, Rsq)
args = rt.RooArgSet(MR, Rsq, nBTaggedJets)
argsWeighted = rt.RooArgSet(MR, Rsq, nBTaggedJets, weight)
importedBackgroundDataSet = rt.RooDataSet("importedBackgroundDataSet", "Combined razor background", argsWeighted, "weight")
for f in range(len(filenames)): 
    if 'SMS' not in filenames[f]: importedBackgroundDataSet.append(rt.RooDataSet(filenames[f], filenames[f], backgroundTrees[f], argsWeighted, "", "weight"))
importedBackgroundDataSet.Print("v")
#rescale weights according to the desired integrated lumi
weightedBackgroundDataSet = rt.RooDataSet("weightedBackgroundDataSet", "Combined razor background", argsWeighted, "weight")
for i in range(0,importedBackgroundDataSet.numEntries()):
   theseArgs = importedBackgroundDataSet.get(i)
   theseArgs.setRealValue("weight", theseArgs.getRealValue("weight")*intLuminosity)
   weightedBackgroundDataSet.add(theseArgs)
weightedBackgroundDataSet.Print("v")
backgroundDataSet3D = convertDataset2Unweighted(weightedBackgroundDataSet, boxName)
bTagString = ""
if bTagBox < 3: bTagString = "nBTaggedJets == "+str(bTagBox)
else: bTagString = "nBTaggedJets >= "+str(bTagBox)
backgroundDataSet = rt.RooDataSet("backgroundDataSet", "Combined razor background", backgroundDataSet3D, args, bTagString)
backgroundDataSet.Print("v")
backgroundDataHist = rt.RooDataHist("backgroundDataHist", "Combined razor background", args2D, backgroundDataSet)
backgroundDataHist.Print("v")

signalDataSet = rt.RooDataSet("signalDataSet", "Razor signal", signalTree, args)
signalDataHist = rt.RooDataHist("signalDataHist", "Razor signal", args2D, signalDataSet)

#define the 2D razor pdf
cutString = "MR > "+str(MRMin)+" && Rsq > "+str(RsqMin)
pars = {}
pars['b'] = rt.RooRealVar("b", "b", 1, 0.0001, 100)
pars['n'] = rt.RooRealVar("n", "n", 3, 0.1, 10)
pars['MR0'] = rt.RooRealVar("MR0", "MR0", -300, -1000, 0.0)
pars['Rsq0'] = rt.RooRealVar("Rsq0", "Rsq0", -.25, -10.0, 0.0)
razorFunction = rt.RooRazor2DTail_SYS("razorFunction", "razorFunction", MR, Rsq, pars['MR0'], pars['Rsq0'], pars['b'], pars['n']);
nTotal = rt.RooRealVar("nTotal", "nTotal", backgroundDataSet.sumEntries(cutString), 0, 10*backgroundDataSet.sumEntries())
extendedRazorFunction = rt.RooExtendPdf("extendedRazorFunction", "Extended razor function", razorFunction, nTotal, "full")

##fit the 2D distribution in each box, plot the distributions and print them to pdfs

###perform fit and plot
        
#define fit region, do fit
fitResult = extendedRazorFunction.fitTo(backgroundDataSet, rt.RooFit.Save(), rt.RooFit.Range("lowMR,lowRsq"), rt.RooFit.Extended(rt.kTRUE))
fitResult.Print("v")

#plot the data and plot the pdf in MR,Rsq
c = rt.TCanvas("c"+boxName, "c"+boxName, 1024, 768)
c.SetLogy()
MRFrame = MR.frame()
backgroundDataSet.plotOn(MRFrame, rt.RooFit.Cut(cutString), rt.RooFit.Binning(MRBinning))  
extendedRazorFunction.plotOn(MRFrame, rt.RooFit.Range("full"), rt.RooFit.LineStyle(2), rt.RooFit.Normalization(nTotal.getVal(), rt.RooAbsReal.NumEvent), rt.RooFit.NormRange("full"), rt.RooFit.ProjectionRange("full"))
RsqFrame = Rsq.frame()
backgroundDataSet.plotOn(RsqFrame, rt.RooFit.Cut(cutString), rt.RooFit.Binning(RsqBinning))
extendedRazorFunction.plotOn(RsqFrame, rt.RooFit.Normalization(nTotal.getVal(), rt.RooAbsReal.NumEvent), rt.RooFit.NormRange("full"), rt.RooFit.ProjectionRange("full"), rt.RooFit.Range("full"), rt.RooFit.LineStyle(2))

#draw and print
MRFrame.SetTitle(boxName+" Box")
MRFrame.SetAxisRange(MRMin, MRMax)
MRFrame.Draw()
printString = outpath+"/RazorExclusionFit_"+boxName+"MR.pdf"
c.Print(printString)

RsqFrame.SetTitle(boxName+" Box")
RsqFrame.SetAxisRange(RsqMin, RsqMax)
RsqFrame.Draw()
printString = outpath+"/RazorExclusionFit_"+boxName+"Rsq.pdf"
c.Print(printString)

###compute the probability of a particular signal, using a Bayesian method

nToys = 10000
pSigmaGivenData = 0 #p(sigma|data) = integral(p(data|sigma,lambda)p(lambda)p(sigma)dlambda)
for toy in range(nToys):
    print("\nRandomizing fit parameters: "+str(toy))
    perturbedArgs = rt.RooArgSet(fitResult.randomizePars())

    B = perturbedArgs.getRealValue("b")
    N = perturbedArgs.getRealValue("n")
    X0 = perturbedArgs.getRealValue("MR0")
    Y0 = perturbedArgs.getRealValue("Rsq0")
    N0 = perturbedArgs.getRealValue("nTotal")

    total_integral = RazorIntegral(B, N, X0, Y0, MRMin, MRMax, RsqMin, RsqMax)
    if total_integral <= 0: #razor function is not well defined in fit region
        print("ERROR: total razor pdf integral = "+str(total_integral))
        continue

    #compute p(data|sigma,lambda)
    logLikelihood = 0
    for [[MRBinMin, MRBinMax],[RsqBinMin, RsqBinMax]] in binCorners: 
        #get estimates for signal and background means in this bin
        meanEvents = RazorIntegral(B, N, X0, Y0, MRBinMin, MRBinMax, RsqBinMin, RsqBinMax)*N0/total_integral
        MR.setVal((MRBinMin+MRBinMax)/2.0)
        Rsq.setVal((RsqBinMin+RsqBinMax)/2.0)
        dataEvents = backgroundDataHist.weight(args2D, 0)
        signalEvents = signalDataHist.weight(args2D, 0)*signalWeight
        print("Data: "+str(dataEvents)+", Bin mean: "+str(meanEvents)+", Signal: "+str(signalEvents))

        #compute the contribution to the log likelihood from this bin (poisson)
        if dataEvents > 0 and meanEvents+signalEvents > 0: 
            logPoisLikelihood = np.log(stats.poisson.pmf(int(dataEvents), meanEvents+signalEvents))
            print(logPoisLikelihood)
            logLikelihood += logPoisLikelihood

    pSigmaGivenData += np.exp(logLikelihood) 

pSigmaGivenData *= pOfSigma
pSigmaGivenData /= nToys

print("Probability of sigma = "+str(sigma)+" given data is "+str(pSigmaGivenData))
