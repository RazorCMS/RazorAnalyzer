import sys, os
import ROOT as rt

def Gamma(a, x):
    return rt.TMath.Gamma(a) * rt.Math.inc_gamma_c(a,x)

def Gfun(x, y, X0, Y0, B, N):
    return Gamma(N,B*N*rt.TMath.Power((x-X0)*(y-Y0),1/N))

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

#create RooDataSets with the MR and Rsq values of events in each tree
MRMin = 400
MRMax = 600
RsqMin = 0.3
RsqMax = 0.6
MRBinWidth = 20
RsqBinWidth = 0.05
MR = rt.RooRealVar("MR", "MR", MRMin, MRMax)
Rsq = rt.RooRealVar("Rsq", "Rsq", RsqMin, RsqMax)
MR.setBins(int((MRMax - MRMin)/MRBinWidth))
Rsq.setBins(int((RsqMax - RsqMin)/RsqBinWidth))
MR.setRange("fullRange", MRMin, MRMax)
Rsq.setRange("fullRange", RsqMin, RsqMax)
MR.setRange("fitRange1", 400, 500)
Rsq.setRange("fitRange1", 0.3, 0.6)
MR.setRange("fitRange2", 400, 600)
Rsq.setRange("fitRange2", 0.3, 0.4)
args = rt.RooArgSet(MR, Rsq)
razorDataSets = dict((box, rt.RooDataSet(box, box, razorBoxes[box], args)) for box in boxNames)
razorDataHists = dict((box, rt.RooDataHist(box, box, args, razorDataSets[box])) for box in boxNames)

#define the 2D razor pdf
b = rt.RooRealVar("b", "b", 1, 0.0001, 100)
n = rt.RooRealVar("n", "n", 3, 0.1, 10)
MR0 = rt.RooRealVar("MR0", "MR0", -300, -1000, 0.0)
Rsq0 = rt.RooRealVar("Rsq0", "Rsq0", -.25, -10.0, 0.0)
razorFunction = rt.RooRazor2DTail_SYS("razorFunction", "razorFunction", MR, Rsq, MR0, Rsq0, b, n);

#fit the 2D distribution in each box, plot the distributions and print them to pdfs
#fit in several different regions and test the reliability of the fit in each region
for box in boxNames:
    #for now only look at hadronic boxes
    if box not in ['MultiJet', '2BJet', '1BJet', '0BJet']: continue
    #if box != 'MultiJet': continue

    #perform fit
    fitResult = razorFunction.fitTo(razorDataSets[box], rt.RooFit.Save(), rt.RooFit.Range("fitRange1"))
    fitResult.Print("v")

    #plot the data and plot the pdf in MR,Rsq
    c = rt.TCanvas("c"+box, "c"+box, 1024, 768)
#    c.SetLogy()
    mRFrame = MR.frame()
    razorDataSets[box].plotOn(mRFrame)  
    razorFunction.plotOn(mRFrame, rt.RooFit.ProjectionRange("fitRange1,fitRange2"), rt.RooFit.NormRange("fitRange1,fitRange2"), rt.RooFit.Range("fullRange"), rt.RooFit.LineStyle(2))
    r2Frame = Rsq.frame()
    razorDataSets[box].plotOn(r2Frame)
    razorFunction.plotOn(r2Frame, rt.RooFit.ProjectionRange("fitRange1,fitRange2"), rt.RooFit.NormRange("fitRange1,fitRange2"), rt.RooFit.Range("fullRange"), rt.RooFit.LineStyle(2))

    #draw and print
    mRFrame.SetTitle(box+" Box")
    mRFrame.Draw()
    printString = outpath+"/RazorInclusiveFits_"+(filename.split("/")[-1]).split(".")[0]+"_"+box+"MR.pdf"
    c.Print(printString)

    r2Frame.SetTitle(box+" Box")
    r2Frame.Draw()
    printString = outpath+"/RazorInclusiveFits_"+(filename.split("/")[-1]).split(".")[0]+"_"+box+"Rsq.pdf"
    c.Print(printString)
