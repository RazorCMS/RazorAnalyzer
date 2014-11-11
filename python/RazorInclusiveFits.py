import sys, os
import string
import numpy as np
import itertools
import ROOT as rt
import rootTools

def Gamma(a, x):
    return rt.TMath.Gamma(a) * rt.Math.inc_gamma_c(a,x)

def Gfun(x, y, X0, Y0, B, N):
    return Gamma(N,B*N*rt.TMath.Power((x-X0)*(y-Y0),1/N))

def RazorIntegral(B, N, X0, Y0, XMin, XMax, YMin, YMax):
    return (N/rt.TMath.Power(B*N,N))*(Gfun(XMin,YMin,X0,Y0,B,N)-Gfun(XMin,YMax,X0,Y0,B,N)-Gfun(XMax,YMin,X0,Y0,B,N)+Gfun(XMax,YMax,X0,Y0,B,N))

##initialization

#check that input file was specified
if len(sys.argv) < 2: 
    print 'Usage: RazorInclusiveFits.py filename.root'
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

print 'Will perform fits on events selected by the razor inclusive analysis'

#load the TTrees from the input file
filename = sys.argv[1]
f = rt.TFile(filename)
boxNames = [
        'MuEle',
        'MuMu',
        'EleEle',
        'MuMultiJet',
        'MuJet',
        'EleMultiJet',
        'EleJet',
        'MultiJet',
        '2BJet',
        '1BJet',
        '0BJet'
        ]
razorBoxes = dict((box, f.Get(box)) for box in boxNames)

##create RooDataSets with the MR and Rsq values of events in each tree

#set bins for MR and Rsq
MRMax = 4000
RsqMax = 1.5
MRBinning = rt.RooBinning(300, MRMax)
RsqBinning = rt.RooBinning(0.15, RsqMax)
MRBinLowEdges = [400, 450, 550, 700, 900, 1200, 1600, 2500]
MRBinHighEdges = [MRBinLowEdges[i] for i in range(1, len(MRBinLowEdges))]
MRBinHighEdges.append(MRMax)
RsqBinLowEdges = [0.2, 0.25, 0.30,0.41,0.52,0.64,0.80]
RsqBinHighEdges = [RsqBinLowEdges[i] for i in range(1, len(RsqBinLowEdges))]
RsqBinHighEdges.append(RsqMax)
for edge in MRBinLowEdges: MRBinning.addBoundary(edge)
for edge in RsqBinLowEdges: RsqBinning.addBoundary(edge)
MR = rt.RooRealVar("MR", "MR", 300, MRMax) 
Rsq = rt.RooRealVar("Rsq", "Rsq", 0.15, RsqMax)
MR.setBinning(MRBinning)
Rsq.setBinning(RsqBinning)

#create several test fit regions
MRMins = [400]
RsqMins = [0.15, 0.2, 0.25]
MRSidebandSizes = [50]
RsqSidebandSizes = [0.05, 0.1, 0.15, 0.2]

#create sideband and extrapolation regions
fitMRMax = 700
fitRsqMax = 0.6 
for MRMin, RsqMin, MRSideband, RsqSideband in itertools.product(MRMins, RsqMins, MRSidebandSizes, RsqSidebandSizes):
    regionName = "m"+str(MRMin)+"r"+string.replace(str(RsqMin), ".", "p")+"ms"+str(MRSideband)+"rs"+string.replace(str(RsqSideband), ".", "p")
    #sidebands (2 rectangles)
    MR.setRange("lowMR"+regionName, MRMin, MRMin+MRSideband)
    Rsq.setRange("lowMR"+regionName, RsqMin, fitRsqMax)
    MR.setRange("lowRsq"+regionName, MRMin, fitMRMax)
    Rsq.setRange("lowRsq"+regionName, RsqMin, RsqMin+RsqSideband)
    #extrapolation region (3 rectangles)
    MR.setRange("extr1"+regionName, fitMRMax, MRMax) #region 1: high MR, low Rsq
    Rsq.setRange("extr1"+regionName, RsqMin, RsqMin+RsqSideband)
    MR.setRange("extr2"+regionName, MRMin, MRMin+MRSideband) #region 2: low MR, high Rsq
    Rsq.setRange("extr2"+regionName, fitRsqMax, RsqMax)
    MR.setRange("extr3"+regionName, MRMin+MRSideband, MRMax) #region 3: high MR, high Rsq
    Rsq.setRange("extr3"+regionName, RsqMin+RsqSideband, RsqMax)
    #full region
    MR.setRange("full"+regionName, MRMin, MRMax)
    Rsq.setRange("full"+regionName, RsqMin, RsqMax)

#create RooDataSets
args = rt.RooArgSet(MR, Rsq)
razorDataSets = dict((box, rt.RooDataSet(box, box, razorBoxes[box], args)) for box in boxNames)
razorDataHists = dict((box, rt.RooDataHist(box, box, args, razorDataSets[box])) for box in boxNames)

#define the 2D razor pdf
pars = {}
pars['b'] = rt.RooRealVar("b", "b", 1, 0.0001, 100)
pars['n'] = rt.RooRealVar("n", "n", 3, 0.1, 10)
pars['MR0'] = rt.RooRealVar("MR0", "MR0", -300, -1000, 0.0)
pars['Rsq0'] = rt.RooRealVar("Rsq0", "Rsq0", -.25, -10.0, 0.0)
razorFunction = rt.RooRazor2DTail_SYS("razorFunction", "razorFunction", MR, Rsq, pars['MR0'], pars['Rsq0'], pars['b'], pars['n']);

##fit the 2D distribution in each box, plot the distributions and print them to pdfs

#fit in several different regions and test the reliability of the fit in each region
for box in boxNames:
    #for now only look at hadronic boxes
    if box not in ['MultiJet', '2BJet', '1BJet', '0BJet']: continue
    #if box != '2BJet': continue

    #test the fit for each sideband region
    for MRMin, RsqMin, MRSideband, RsqSideband in itertools.product(MRMins, RsqMins, MRSidebandSizes, RsqSidebandSizes):
        #MR in the multijet box peaks around 350, so don't try fitting if MRMin < 400
        if box == 'MultiJet' and MRMin < 399: continue
        print("MRMin = "+str(MRMin)+", RsqMin = "+str(RsqMin))

        ###perform fit and plot
        
        #reset fit parameters
        pars['b'].setVal(1.0)
        pars['n'].setVal(3.0)
        pars['MR0'].setVal(-300)
        pars['Rsq0'].setVal(-.25)

        #define fit region, do fit
        regionName = "m"+str(MRMin)+"r"+string.replace(str(RsqMin), ".", "p")+"ms"+str(MRSideband)+"rs"+string.replace(str(RsqSideband), ".", "p")
        fitResult = razorFunction.fitTo(razorDataSets[box], rt.RooFit.Save(), rt.RooFit.Range("lowMR"+regionName+",lowRsq"+regionName))
        fitResult.Print("v")

        #plot the data and plot the pdf in MR,Rsq
        c = rt.TCanvas("c"+box, "c"+box, 1024, 768)
        c.SetLogy()
        MRFrame = MR.frame()
        cutString = "MR > "+str(MRMin)+" && Rsq > "+str(RsqMin)
        razorDataSets[box].plotOn(MRFrame, rt.RooFit.Cut(cutString), rt.RooFit.Binning(MRBinning))  
        razorFunction.plotOn(MRFrame, rt.RooFit.Normalization(razorDataSets[box].sumEntries(cutString), rt.RooAbsReal.NumEvent), rt.RooFit.NormRange("full"+regionName), rt.RooFit.ProjectionRange("full"+regionName), rt.RooFit.Range("full"+regionName), rt.RooFit.LineStyle(2))
        RsqFrame = Rsq.frame()
        razorDataSets[box].plotOn(RsqFrame, rt.RooFit.Cut(cutString), rt.RooFit.Binning(RsqBinning))
        razorFunction.plotOn(RsqFrame, rt.RooFit.Normalization(razorDataSets[box].sumEntries(cutString), rt.RooAbsReal.NumEvent), rt.RooFit.NormRange("full"+regionName), rt.RooFit.ProjectionRange("full"+regionName), rt.RooFit.Range("full"+regionName), rt.RooFit.LineStyle(2))

        #draw and print
        MRFrame.SetTitle(box+" Box, "+regionName)
        MRFrame.SetAxisRange(MRMin, MRMax)
        MRFrame.Draw()
        printString = outpath+"/RazorInclusiveFits_"+(filename.split("/")[-1]).split(".")[0]+"_"+box+"MR"+regionName+".pdf"
        c.Print(printString)

        RsqFrame.SetTitle(box+" Box, "+regionName)
        RsqFrame.SetAxisRange(RsqMin, RsqMax)
        RsqFrame.Draw()
        printString = outpath+"/RazorInclusiveFits_"+(filename.split("/")[-1]).split(".")[0]+"_"+box+"Rsq"+regionName+".pdf"
        c.Print(printString)

        ###validation of the fit

        #perform a series of toy experiments: 
        #vary fit function according to the fit covariance matrix and compute the variance of each fit yield
        nToys = 100
        fitArgs = fitResult.floatParsFinal()
        bestFitArgs = 
        binCorners = [[[MRBinMin,MRBinMax], [RsqBinMin, RsqBinMax]] for [(MRBinMin,MRBinMax), (RsqBinMin,RsqBinMax)] in itertools.product(zip(MRBinLowEdges,MRBinHighEdges), zip(RsqBinLowEdges,RsqBinHighEdges)) if MRBinMin+1 >= fitMRMax or RsqBinMin+0.01 >= fitRsqMax or (MRBinMin+1 >= MRMin+MRSideband and RsqBinMin+0.01 >= RsqMin+RsqSideband)] #list of all bins in the extrapolation region
        binEvents = []
        NTOT = razorDataSets[box].sumEntries(cutString) #total number of events in fit range + extrapolation region
        for toy in range(nToys):
            #first time through: fit parameters unperturbed

            #get parameter values
            B = pars['b'].getVal()
            N = pars['n'].getVal()
            X0 = pars['MR0'].getVal()
            Y0 = pars['Rsq0'].getVal()

            #check if the fit function is sane by integrating it over the whole range
            total_integral = RazorIntegral(B, N, X0, Y0, MRMin, MRMax, RsqMin, RsqMax)
            if total_integral <= 0: #razor function is not well defined
                print("\nERROR: total razor pdf integral ="+str(total_integral))
                continue

            #in each bin in the extrapolation region, integrate the function to get the expected number of events
            integrals = [RazorIntegral(B, N, X0, Y0, MRBinMin, MRBinMax, RsqBinMin, RsqBinMax) for [[MRBinMin,MRBinMax],[RsqBinMin,RsqBinMax]] in binCorners]
            theseBinEvents = [NTOT*integral/total_integral for integral in integrals]
            binEvents.append(theseBinEvents)

            #perturb fit (assumes a gaussian likelihood!) and set fit params equal to the perturbed values
            perturbedArgs = fitResult.randomizePars()
            for p in rootTools.RootIterator.RootIterator(perturbedArgs): pars[p.GetName()].setVal(p.getVal())

        if len(binEvents) == 0: continue
        #get the number of events and variance in each bin
        binContentsFit = binEvents[0]
        sigmaFit = [np.var(np.array([binEvents[toy][bin] for toy in range(1,len(binEvents))])) for bin in range(len(binCorners))]

        #get the number of data points in each bin
        binWeights = []
        for [[MRBinMin,MRBinMax],[RsqBinMin,RsqBinMax]] in binCorners:
            MR.setVal((MRBinMin+MRBinMax)/2.0) #input coordinates of a point in the bin
            Rsq.setVal((RsqBinMin+RsqBinMax)/2.0)
            binWeight = razorDataHists[box].weight(args, 0) #0 -> no interpolation
            binWeights.append(binWeight)

        #compute the pull in each bin
        pulls = rt.TH1F("pulls"+regionName, "Pull distribution, "+box+" Box, "+regionName+"; (data - fit)/#sigma_{fit}", 20, -4, 4)
        for bin in range(len(binCorners)): pulls.Fill((binWeights[bin] - binContentsFit[bin])/sigmaFit[bin])

        #plot the distribution of pulls, fit with a gaussian, and print to pdf

        cPull = rt.TCanvas("cPull"+regionName, "cPull"+regionName, 800, 800)
        pulls.Fit("gaus")
        pulls.Draw()
        cPull.Print(outpath+"/razorInclusivePulls_"+(filename.split("/")[-1]).split(".")[0]+"_"+box+regionName+".pdf")

        ###validation 2: mixed-sample test

        n_k = 10

        #our dataset, restricted to extrapolation region
        extrRegionString = "MR > "+fitMRMax+" || Rsq > "+fitRsqMax+" || (MR > "+MRMin+MRSideband+" && Rsq > "+RsqMin+RsqSideband+")"
        razorDataExtr = rt.RooDataSet("razorDataExtr", "razorDataExtr", razorDataSets[box], args, extrRegionString)
        nExtr = razorDataExtr.sumEntries()

        #we need to create a sample of size genSampleSize whose events are all in the signal region.
        for p in rootTools.RootIterator.RootIterator(fitArgs): pars[p.GetName()].setVal(p.getVal())
        genSpec = razorFunction.prepareMultiGen(args, rt.RooFit.NumEvent(1)) #prepare pdf for generation
        genSampleSize = 10*nExtr #size of the toy dataset used to compute the test statistic
        genSample = rt.RooDataSet("genSample", "genSample", args) #empty dataset
        while genSample.sumEntries() < genSampleSize:
            genDataPoint = razorFunction.generate(genSpec)
            genArgs = genDataPoint.get(0)
            if genArgs.getRealValue("MR") > fitMRMax or genArgs.getRealValue("Rsq") > fitRsqMax or (genArgs.getRealValue("MR") > MRMin + MRSideband and genArgs.getRealValue("Rsq") > RsqMin + RsqSideband): genSample.append(genDataPoint)

        #use the MLMixedSample object to get the test statistic value
        mlMixed = rt.MLMixedSample(n_k);
        testStat = mlMixed.testStatistic(razorDataExtr, genSample, args)

        #compute the p-value associated with the test statistic
        n_a = nExtr
        n_b = genSampleSize
        mu_T = (n_a*(n_a-1) + n_b*(n_b-1))*1.0/((n_a+n_b)*(n_a+n_b-1))
        #limiting form of the variance; should be accurate for large sample sizes
        var_T = 1.0/((n_a+n_b)*n_k)*(n_a*n_b*1.0/((n_a+n_b)*(n_a+n_b))+4*n_a*n_a*n_b*n_b*1.0/((n_a+n_b)*(n_a+n_b)*(n_a+n_b)*(n_a+n_b)))
        pValue = 0.5*(1 - erf((testStat - mu_T)/(var_T*sqrt(2))))
        print("The mixed-sample test statistic is "+str(testStat)+".  The mean is "+str(mu_T)+" and the variance is "+str(var_T)+".  This measurement thus yields a p-value of "+str(pValue))
