import sys, os
import string
import itertools
import ROOT as rt

def Gamma(a, x):
    return rt.TMath.Gamma(a) * rt.Math.inc_gamma_c(a,x)

def Gfun(x, y, X0, Y0, B, N):
    return Gamma(N,B*N*rt.TMath.Power((x-X0)*(y-Y0),1/N))

##initialization

#check that input file was specified
if len(sys.argv) < 2: 
    print 'Usage: RazorInclusiveFits.py filename.root'
    exit()

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
razorTree = f.Get("RazorInclusive")

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

##create RooDataSets with the MR and Rsq values of events in each tree

#set number of bins for MR and Rsq
MRBinWidth = 20
RsqBinWidth = 0.05
MRMax = 1000
RsqMax = 0.8
MR = rt.RooRealVar("MR", "MR", 300, MRMax) 
Rsq = rt.RooRealVar("Rsq", "Rsq", 0.15, RsqMax)
box = rt.RooRealVar("box", "box", 0, 100)
MR.setBins(int(MRMax/MRBinWidth))
Rsq.setBins(int(RsqMax/RsqBinWidth))

#create several test fit regions
MRMins = [300, 350, 400, 450]
RsqMins = [0.15, 0.2, 0.25]
MRSidebandSizes = [100, 150, 200]
RsqSidebandSizes = [0.1, 0.15, 0.2]

#create sideband and extrapolation regions
fitMRMax = 600
fitRsqMax = 0.6 
for MRMin, RsqMin, MRSideband, RsqSideband in itertools.product(MRMins, RsqMins, MRSidebandSizes, RsqSidebandSizes):
    regionName = "m"+str(MRMin)+"r"+string.replace(str(RsqMin), ".", "p")+"ms"+str(MRSideband)+"rs"+string.replace(str(RsqSideband), ".", "p")
    #sidebands (2 rectangles)
    MR.setRange("lowMR"+regionName, MRMin, MRMin+MRSideband)
    Rsq.setRange("lowMR"+regionName, RsqMin, 0.6)
    MR.setRange("lowRsq"+regionName, MRMin, 600)
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

args = rt.RooArgSet(MR, Rsq, box)
razorDataSets = dict((box, rt.RooDataSet(box, box, razorTree, args, "box == "+str(boxNames.index(box)))) for box in boxNames)
razorDataHists = dict((box, rt.RooDataHist(box, box, args, razorDataSets[box])) for box in boxNames)

#define the 2D razor pdf
b = rt.RooRealVar("b", "b", 1, 0.0001, 100)
n = rt.RooRealVar("n", "n", 3, 0.1, 10)
MR0 = rt.RooRealVar("MR0", "MR0", -300, -1000, 0.0)
Rsq0 = rt.RooRealVar("Rsq0", "Rsq0", -.25, -10.0, 0.0)
razorFunction = rt.RooRazor2DTail_SYS("razorFunction", "razorFunction", MR, Rsq, MR0, Rsq0, b, n);

##fit the 2D distribution in each box, plot the distributions and print them to pdfs

#fit in several different regions and test the reliability of the fit in each region
for box in boxNames:
    #for now only look at hadronic boxes
    #if box not in ['MultiJet', '2BJet', '1BJet', '0BJet']: continue
    if box != '2BJet': continue

    #test the fit for each sideband region
    for MRMin, RsqMin, MRSideband, RsqSideband in itertools.product(MRMins, RsqMins, MRSidebandSizes, RsqSidebandSizes):
        #MR in the multijet box peaks around 350, so don't try fitting if MRMin < 400
        if box == 'MultiJet' and MRMin < 399: continue

        #perform fit
        regionName = "m"+str(MRMin)+"r"+string.replace(str(RsqMin), ".", "p")+"ms"+str(MRSideband)+"rs"+string.replace(str(RsqSideband), ".", "p")
        fitResult = razorFunction.fitTo(razorDataSets[box], rt.RooFit.Save(), rt.RooFit.Range("lowMR"+regionName+",lowRsq"+regionName))
        fitResult.Print("v")

        #plot the data and plot the pdf in MR,Rsq
        c = rt.TCanvas("c"+box, "c"+box, 1024, 768)
        #c.SetLogy()
        mRFrame = MR.frame()
        #razorDataSets[box].plotOn(mRFrame)  
        razorDataSets[box].plotOn(mRFrame, rt.RooFit.Cut("Rsq > "+str(RsqMin)))  
        razorFunction.plotOn(mRFrame, rt.RooFit.ProjectionRange("lowMR"+regionName+",lowRsq"+regionName), rt.RooFit.NormRange("lowMR"+regionName+",lowRsq"+regionName), rt.RooFit.Range("full"+regionName), rt.RooFit.LineStyle(2))
        r2Frame = Rsq.frame()
        #razorDataSets[box].plotOn(r2Frame)
        razorDataSets[box].plotOn(r2Frame, rt.RooFit.Cut("MR > "+str(MRMin)))
        razorFunction.plotOn(r2Frame, rt.RooFit.ProjectionRange("lowMR"+regionName+",lowRsq"+regionName), rt.RooFit.NormRange("lowMR"+regionName+",lowRsq"+regionName), rt.RooFit.Range("full"+regionName), rt.RooFit.LineStyle(2))

        #draw and print
        mRFrame.SetTitle(box+" Box, "+regionName)
        mRFrame.SetAxisRange(MRMin, MRMax)
        mRFrame.Draw()
        printString = outpath+"/RazorInclusiveFits_"+(filename.split("/")[-1]).split(".")[0]+"_"+box+"MR"+regionName+".pdf"
        c.Print(printString)

        r2Frame.SetTitle(box+" Box, "+regionName)
        r2Frame.SetAxisRange(RsqMin, RsqMax)
        r2Frame.Draw()
        printString = outpath+"/RazorInclusiveFits_"+(filename.split("/")[-1]).split(".")[0]+"_"+box+"Rsq"+regionName+".pdf"
        c.Print(printString)
