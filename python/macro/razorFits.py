import os
import numpy as np
import ROOT as rt 
from array import array

# Local imports
from razorAnalysis import Analysis
import macro
from framework import Config
from DustinTuple2RooDataSet import convertTree2Dataset, initializeWorkspace
from WriteDataCard import convertDataset2TH1
from WriteDataCard import initializeWorkspace as initializeFitWorkspace
from PlotFit import setStyle, convertSideband, get3DHistoFrom1D, print1DProj, get1DHistoFrom2D, print2DResiduals
from BinnedFit import binnedFit

def unweightHist(hist):
    """Throws a Poisson toy from each bin of the histogram.
        Returns the resulting unweighted histogram"""
    newHist = hist.Clone()
    for ibin in range(1, newHist.GetNbinsX()+1):
        newHist.SetBinContent( ibin, np.random.poisson(
            hist.GetBinContent(ibin)) )
    return newHist

def make3DHistProjections(hist3D):
    """Makes XY, X, and Y projections of hist3D.
        Returns a list containing the three projections."""
    x = hist3D.Project3D("xe")
    y = hist3D.Project3D("ye")
    xy = hist3D.Project3D("yxe")
    return [x,y,xy]

def getStatNSigmaHist(fitHist, datHist, name):
    """Return a histogram with (data-fit)/sqrt(fit) in each bin.
        fitHist: 2D histogram with fit predictions
        datHist: 2D histogram with data yields"""
    statNSigma = fitHist.Clone(name)
    for i in range(1,statNSigma.GetNbinsX()+1):
        for j in range(1,statNSigma.GetNbinsY()+1):
            fit = fitHist.GetBinContent(i,j)
            dat = datHist.GetBinContent(i,j)
            if fit > 0.0:
                statNSigma.SetBinContent(i,j,(dat-fit)/rt.TMath.Sqrt(fit))                        
            else:
                print "ERROR FIT = 0, SKIPPING BIN"
    return statNSigma

def getBtagLabel(binning):
    btagLabel = ""
    if binning[-1] == binning[0]+1 and binning[-1]==4:
        btagLabel = "#geq %i b-tag" % binning[0]
    elif binning[-1] == binning[0]+1:
        btagLabel = "%i b-tag" % binning[0]
    elif binning[-1]==4:
        btagLabel = "#geq %i b-tag" % binning[0]
    elif len(binning)==2 and binning[0]==1 and binning[-1]==4:
        btagLabel = "#geq %i b-tag" % binning[0]        
    else:
        btagLabel = "%i-%i b-tag" % (binning[0],binning[-2])
    return btagLabel

class FitError(Exception):
    """Exception for fit-related errors"""
    pass

class FitInstance(object):
    """Helper class for doing razor fit bookkeeping.
        attributes:
          analysis: Analysis object storing cuts and filenames
          weights: dict of the form { process:weight }
          isData: bool
          workspace: RooWorkspace to hold dataset, fit result, fit function, etc
          filename: name of ROOT file used to hold workspace
          config: Config object defining fit function parameters and binning
          fitRegion: either 'Full' or 'LowMR,LowRsq'
          sideband: string containing sideband bin information
          x,y,z: arrays with bin boundaries
    """

    def __init__(self, box, tag="Razor2016_MoriondRereco", 
            weights=None, isData=True, configFile="config/run2_2016.config",
            full=False):
        self.sideband = None
        self.weights = weights
        self.isData = isData
        if full:
            self.fitRegion = 'Full'
        else:
            self.fitRegion = 'LowMR,LowRsq'
        self.analysis = Analysis(box, tag=tag)
        self.filename = self.getFilename()
        self.config = Config.Config(configFile)
        self.x = array('d', self.config.getBinning(box)[0]) # MR binning
        self.y = array('d', self.config.getBinning(box)[1]) # Rsq binning
        self.z = array('d', self.config.getBinning(box)[2]) # nBtag binning
        self.workspace = rt.RooWorkspace("w"+box)
        initializeWorkspace(self.workspace, self.config, box)
        setStyle()

    def addToWorkspace(self, obj, tobject=False):
        """Accesses the appropriate import method of RooWorkspace.
            Set tobject=True to import as a generic TObject."""
        # A blank RooCmdArg is needed so that the correct overloaded 
        # import() method is used.
        method = getattr(self.workspace,'import')
        if not tobject:
            method(obj, rt.RooCmdArg())
        else:
            method(obj)

    def getFilename(self):
        """Builds the name of the file to use for input/output"""
        strs = ["RazorFitInstance",self.analysis.tag,self.analysis.region]
        if self.fitRegion == 'Full':
            strs.append('Full')
        if not self.isData:
            strs.append("SMCocktail")
            if self.weights is not None:
                for proc,weight in self.weights.iteritems():
                    strs.append(("%s_%.2f"%(proc,weight)).replace(".","p"))
        return "_".join(strs)+".root"

    def getTrees(self):
        """Opens the ntuple files and loads the tree from each"""
        if self.isData:
            filenames = {"Data":self.analysis.filenames["Data"]}
        else:
            filenames = {proc:self.analysis.filenames[proc] for proc in 
                    self.analysis.filenames if proc != "Data" and proc != "QCD"}
        files = { proc:rt.TFile.Open(filenames[proc]) for proc in filenames }
        trees = macro.makeTreeDict(files, "RazorInclusive")
        return trees

    def initDataset(self):
        """Loads the ntuple files and creates the dataset for fitting"""
        trees = self.getTrees()
        datasets = []
        for proc,tree in trees.iteritems():
            print proc+":"
            scale = 1.0
            if proc in self.weights:
                scale = self.weights[proc]
                print "Scaling yields by %.2f"%(scale)
            datasets.append( convertTree2Dataset(tree, self.config,
                self.analysis.region, self.workspace, 
                useWeight=not self.isData, 
                globalScaleFactor=self.analysis.lumi*scale,
                isData=self.isData, tag=self.analysis.tag) )
        combinedData = datasets[0].Clone("RMRTree")
        for i in range(1,len(datasets)):
            combinedData.append(datasets[i])
        self.addToWorkspace(combinedData)

    def initBinnedDataset(self):
        """Initializes the binned razor pdf and dataset"""
        initializeFitWorkspace(self.workspace, self.config, self.analysis.region)
        nBins = (len(self.x)-1)*(len(self.y)-1)*(len(self.z)-1)
        th1x = self.workspace.var('th1x')
        th1x.setBins(nBins)
        hist = convertDataset2TH1(self.workspace.data('RMRTree'), self.config, 
                self.analysis.region, self.workspace)
        if not self.isData:
            hist = unweightHist(hist)
        dataHist = rt.RooDataHist("data_obs","data_obs",rt.RooArgList(th1x), 
                rt.RooFit.Import(hist))
        dataHist.Print('v')
        self.addToWorkspace(dataHist)

    def loadWorkspace(self, filename=None):
        """Loads the workspace from file.  If no filename is given,
            the name returned by get_filename() will be used"""
        if filename is None:
            filename = self.filename
        f = rt.TFile.Open(filename)
        if not f:
            raise FitError("File %s could not be opened"%(filename))
        wName = "w"+self.analysis.region
        self.workspace = f.Get(wName)
        if not self.workspace:
            raise FitError("Workspace %s not found in file %s"%(
                wName, filename))
        print "Loaded workspace %s from file %s"%(wName,filename)

    def writeWorkspace(self, filename=None):
        """Saves the workspace to a file.  If no filename is given,
            the name returned by get_filename() will be used"""
        if filename is None:
            filename = self.filename
        f = rt.TFile(filename, "RECREATE")
        self.workspace.Write()
        print "Wrote workspace to file %s"%(filename)
        f.Close()

    def fit(self):
        """Fits the razor pdf to the data"""
        extRazorPdf = self.workspace.pdf('extRazorPdf')
        datahist = self.workspace.data('data_obs')
        self.sideband = convertSideband(self.fitRegion, self.workspace, 
                self.x, self.y, self.z)
        result = binnedFit(extRazorPdf, datahist, self.sideband)
        result.Print('v')
        self.addToWorkspace(result, tobject=True)

    def get3DFitHistos(self):
        """Returns a pair of TH3s: (data, fit prediction)"""
        extRazorPdf = self.workspace.pdf('extRazorPdf')
        th1x = self.workspace.var('th1x')
        dataHist = self.workspace.data("data_obs")
        asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),
                rt.RooFit.Name('central'),rt.RooFit.Asimov())
        opt = [rt.RooFit.CutRange(myRange) 
                for myRange in self.sideband.split(',')]
        asimov_reduce = asimov.reduce(opt[0])
        dataHist_reduce = dataHist.reduce(opt[0])
        for iOpt in range(1,len(opt)):
            asimov_reduce.add(asimov.reduce(opt[iOpt]))
            dataHist_reduce.add(dataHist.reduce(opt[iOpt]))
        h_th1x = asimov_reduce.createHistogram('h_th1x',th1x)
        h_data_th1x = dataHist_reduce.createHistogram('h_data_th1x',th1x)
        h_data_nBtagRsqMR = get3DHistoFrom1D(h_data_th1x,self.x,self.y,self.z,"h_data_nBtagRsqMR")
        h_nBtagRsqMR = get3DHistoFrom1D(h_th1x,self.x,self.y,self.z,"h_nBtagRsqMR")
        return h_data_nBtagRsqMR, h_nBtagRsqMR

    def getSidebandMax(self):
        sideband = None
        if self.fitRegion=="LowMR,LowRsq":
            mrSide = self.workspace.var('MR').getMax('LowMR')
            rsqSide = self.workspace.var('Rsq').getMax('LowRsq')
            sideband = [mrSide, rsqSide]
        return sideband

    def plot(self, filename=None):
        """Plots the fit results (before toys)"""
        if filename is None:
            filename = self.filename.replace('.root','_Plots.root')
        f = rt.TFile(filename, 'UPDATE')
        dirName = "BeforeToys"
        tdirectory = f.GetDirectory(dirName)
        if tdirectory==None:
            f.mkdir(dirName)
            tdirectory = f.GetDirectory(dirName)
        c = rt.TCanvas('c','c',500,400)
        rt.SetOwnership(c, False)
        rt.TH1D.SetDefaultSumw2()
        rt.TH2D.SetDefaultSumw2()
        rt.TH3D.SetDefaultSumw2()
        if self.sideband is None:
            self.sideband = convertSideband(self.fitRegion, self.workspace, 
                    self.x, self.y, self.z)

        h_data_nBtagRsqMR, h_nBtagRsqMR = self.get3DFitHistos()
        data_projections = make3DHistProjections(h_data_nBtagRsqMR)
        fit_projections = make3DHistProjections(h_nBtagRsqMR)
        for h in data_projections+fit_projections:
            tdirectory.cd()
            h.Write()
        if len(self.z)>1:
            h_MR_components = []
            h_Rsq_components = []
            h_RsqMR_components = []
            h_data_RsqMR_components = []
            h_sig_RsqMR_components = []
            h_th1x_components = []
            h_data_th1x_components = []
            h_sig_th1x_components = []
            h_labels = []        
            h_colors = [rt.kOrange,rt.kViolet,rt.kRed,rt.kGreen,rt.kGray+2]
            for k in range(1,len(self.z)):
                h_MR_components.append(h_nBtagRsqMR.ProjectionX(
                    "MR_%ibtag"%self.z[k-1],0,-1,k,k,""))
                h_Rsq_components.append(h_nBtagRsqMR.ProjectionY(
                    "Rsq_%ibtag"%self.z[k-1],0,-1,k,k,""))
                h_nBtagRsqMR.GetZaxis().SetRange(k,k)
                h_RsqMR_components.append(h_nBtagRsqMR.Project3D(
                    "%ibtag_yx"%self.z[k-1]))
                h_data_nBtagRsqMR.GetZaxis().SetRange(k,k)
                h_data_RsqMR_components.append(h_data_nBtagRsqMR.Project3D(
                    "%ibtag_yx"%self.z[k-1]))
                h_th1x_components.append(get1DHistoFrom2D(
                    h_RsqMR_components[-1],self.x,self.y,
                    'h_th1x_%ibtag'%(self.z[k-1])))
                h_data_th1x_components.append(get1DHistoFrom2D(
                    h_data_RsqMR_components[-1],self.x,self.y,
                    'h_th1x_data_%ibtag'%(self.z[k-1])))
                if self.z[k-1]==3 and self.z[-1]==4:
                    h_labels.append("#geq %i b-tag" % self.z[k-1] )
                elif self.z[k-1]==1 and self.z[-1]==4 and len(self.z)==2:
                    h_labels.append("#geq %i b-tag" % self.z[k-1] )
                else:            
                    h_labels.append("%i b-tag" % self.z[k-1] )

        h_RsqMR_statnsigma_components = []
        if len(self.z)>1:
            for k in range(1,len(self.z)):
                h_RsqMR_statnsigma_components.append( 
                        getStatNSigmaHist(h_RsqMR_components[k-1], 
                            h_data_RsqMR_components[k-1],
                            "h_RsqMR_statnsigma_%ibtag"%self.z[k-1]))

        btagLabel = getBtagLabel(self.z)
        lumiLabel = "%.0f fb^{-1} (13 TeV)" % (self.analysis.lumi/1000.)
        boxLabel = "razor %s %s %s Fit" % (self.analysis.region,
                btagLabel,self.fitRegion.replace('LowMR,LowRsq','Sideband'))
        plotLabel = ""
        sidebandFit = self.getSidebandMax()
        for h in h_RsqMR_components:
            tdirectory.cd()
            h.Write()
        plotDir = "Plots/%s/Fits"%(self.analysis.tag)
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)
        print1DProj(c,tdirectory,fit_projections[0],data_projections[0],
                plotDir+"/h_MR_%s.pdf"%self.analysis.region,"M_{R} [GeV]",
                "Events",lumiLabel,boxLabel,plotLabel,self.isData,False,None,
                None,h_MR_components,h_colors,h_labels)
        print1DProj(c,tdirectory,fit_projections[1],data_projections[1],
                plotDir+"/h_Rsq_%s.pdf"%self.analysis.region,"R^{2}",
                "Events",lumiLabel,boxLabel,plotLabel,self.isData,False,None,
                None,h_Rsq_components,h_colors,h_labels)
        if len(self.z)>2:
            for k in range(0,len(self.z)-1):
                newBoxLabel = "razor %s %s %s Fit"%(self.analysis.region,
                    h_labels[k],self.fitRegion.replace('LowMR,LowRsq','Sideband'))
                print1DProj(c,tdirectory,h_th1x_components[k],
                    h_data_th1x_components[k], 
                    plotDir+"/h_th1x_%ibtag_%s.pdf"%(self.z[k],
                        self.analysis.region),"Bin Number",
                        "Events",lumiLabel,newBoxLabel,plotLabel,self.isData,
                        False, None)
                print2DResiduals(c,tdirectory,h_RsqMR_statnsigma_components[k],
                    plotDir+"/h_RsqMR_statnsigma_log_%ibtag_%s.pdf"%(self.z[k],
                    self.analysis.region), "M_{R} [GeV]", "R^{2}", 
                    "Stat. n#sigma (Data - Fit)/sqrt(Fit)",
                    lumiLabel,newBoxLabel,plotLabel,self.x,self.y,self.isData,
                    sidebandFit,False,None)
        f.Close()

    def doFitSequence(self, load=False, doFit=True, plot=True):
        """Performs all steps needed to build and fit the dataset"""
        if load:
            self.loadWorkspace()
        else:
            self.initDataset()
            self.initBinnedDataset()
        if doFit:
            self.fit()
        if plot:
            self.plot()
        self.writeWorkspace()
