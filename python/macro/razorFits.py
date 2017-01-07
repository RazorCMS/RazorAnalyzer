import os
import numpy as np
import ROOT as rt 
from array import array

# Local imports
import macro
from razorAnalysis import Analysis
from framework import Config
from DustinTuple2RooDataSet import convertTree2Dataset, initializeWorkspace
from BinnedFit import binnedFit
from WriteDataCard import convertDataset2TH1
from WriteDataCard import initializeWorkspace as initializeFitWorkspace
from RunToys import runToys
from PlotFit import setStyle, convertSideband, get3DHistoFrom1D, print1DProj, print1DProjNs, get1DHistoFrom2D, print2DResiduals, getErrors1D, getErrors2D, getErrors3D, getNsigma2D

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
          dirname: name of directory for output files
          toysFile: path to file containing stat+sys toys
          sysFile: path to file containing sys-only toys
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
        self.dirname = os.path.dirname(self.filename)
        self.toysFile = self.getToysFilename(noStat=False)
        self.sysFile = self.getToysFilename(noStat=True)
        if not os.path.isdir(self.dirname):
            os.makedirs(self.dirname)
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
            method(obj, True)

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
        filename = "_".join(strs)+".root"
        return self.getDirname(filename)+"/"+filename

    def getDirname(self, filename):
        """Builds the name of the output directory"""
        dirname = "Plots/%s/Fits_%s"%(self.analysis.tag, 
                filename.replace(".root","").replace(
                "RazorFitInstance_%s_"%self.analysis.tag,""))
        return dirname

    def getToysFilename(self, noStat=False):
        """Gets the appropriate toys filename"""
        if noStat:
            unc = "noStat"
        else:
            unc = "varyN"
        fname = "toys_Bayes_%s_%s.root"%(unc, self.analysis.region)
        return self.dirname+"/"+fname

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

    def setDefaultFitParams(self):
        """Gets the default fit parameters from the config and applies
            them to the variables in the current workspace"""
        for param in self.config.getVariables(self.analysis.region, 
                "combine_parameters"):
            name,rest = param.replace(']','').split('[')
            val = float(rest.split(',')[0])
            if not ("Cut" in name or "Ntot" in name):
                print "Resetting fit parameter %s to %.2f"%(name,val)
                self.workspace.var(name).setVal(val)
        # Set Ntot parameters according to number of events in dataset
        for k in range(0, len(self.z)-1):
            entries = self.workspace.data("RMRTree").sumEntries(
                    "nBtag>=%i && nBtag<%i"%(self.z[k], self.z[k+1]))
            par = "Ntot_TTj%db_MultiJet"%self.z[k]
            print "Resetting fit parameter %s to %d"%(par, entries)
            self.workspace.var(par).setVal(entries)

    def fit(self):
        """Fits the razor pdf to the data"""
        self.setDefaultFitParams()
        extRazorPdf = self.workspace.pdf('extRazorPdf')
        datahist = self.workspace.data('data_obs')
        self.sideband = convertSideband(self.fitRegion, self.workspace, 
                self.x, self.y, self.z)
        result = binnedFit(extRazorPdf, datahist, self.sideband)
        result.Print('v')
        self.addToWorkspace(result, tobject=True)

    def plotCorrelationMatrix(self):
        """Plots the correlation matrix of the fit parameters"""
        c = rt.TCanvas('c','c',400,300)
        self.workspace.obj("nll_extRazorPdf_data_obs").correlationHist(
                ).Draw("colz")
        c.Print(self.dirname+"/correlationHist.pdf")

    def get3DFitHistos(self, sideband):
        """Returns a pair of TH3s: (data, fit prediction)"""
        extRazorPdf = self.workspace.pdf('extRazorPdf')
        th1x = self.workspace.var('th1x')
        dataHist = self.workspace.data("data_obs")
        asimov = extRazorPdf.generateBinned(rt.RooArgSet(th1x),
                rt.RooFit.Name('central'),rt.RooFit.Asimov())
        opt = [rt.RooFit.CutRange(myRange) 
                for myRange in sideband.split(',')]
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
        """Returns a list containing the upper edge values of the
            MR and Rsq sidebands"""
        sideband = None
        if self.fitRegion=="LowMR,LowRsq":
            mrSide = self.workspace.var('MR').getMax('LowMR')
            rsqSide = self.workspace.var('Rsq').getMax('LowRsq')
            sideband = [mrSide, rsqSide]
        return sideband

    def runToys(self, sysPlusStat=True, numToys=1000):
        """Call the toy-generation machinery"""
        # This is a dummy class to hold options for runToys()
        class Options(object):
            pass
        opts = Options()
        opts.r = -1
        opts.box = self.analysis.region
        opts.fitRegion = self.fitRegion
        opts.varyN = sysPlusStat
        opts.noStat = (not sysPlusStat)
        opts.noSys = False
        opts.freq = False
        opts.oneSigma = False
        opts.outDir = self.dirname
        opts.nToys = numToys
        toysFile = runToys(self.workspace, opts, self.config, seed=-1)
        print "Uncertainty information stored in %s"%toysFile
        return toysFile

    def getToyTree(self, inputFile=None):
        """Get tree from file"""
        toyTree = None
        if inputFile is not None:
            toyFiles = inputFile.split(',')
            toyTree = rt.TChain("myTree")
            for toyFile in toyFiles:
                toyTree.Add(toyFile)
            toyTree.SetName("toyTree")
        return toyTree

    def blindSideband(self, hist, plotRegion):
        w = self.workspace
        if plotRegion!='Full':            
            for i in range(1,hist.GetNbinsX()+1):
                for j in range(1,hist.GetNbinsY()+1):
                    w.var('MR').setVal(hist.GetXaxis().GetBinCenter(i))
                    w.var('Rsq').setVal(hist.GetYaxis().GetBinCenter(j))
                    w.var('nBtag').setVal(self.z[k-1])
                    inSideband = 0
                    for fitname in plotRegion.split(','):
                        inSideband += ( w.var('MR').inRange(fitname) 
                                * w.var('Rsq').inRange(fitname) 
                                * w.var('nBtag').inRange(fitname) )
                    if not inSideband:
                        hist.SetBinContent(i,j,0)

    def getPlottingOptions(self):
        """Wraps appropriate plotting options in an object and returns it"""
        class Options(object):
            pass
        options = Options()
        options.fitRegion = self.fitRegion
        options.noStat = True
        options.printErrors = True
        options.outDir = self.dirname
        return options

    def plot(self, filename=None, unblind=False, 
            toysFile=None, sysFile=None):
        """Plots the fit results"""
        options = self.getPlottingOptions()
        if unblind:
            plotRegion = 'Full'
        else:
            plotRegion = self.fitRegion
        if filename is None:
            filename = self.filename.replace('.root','_Plots.root')
        f = rt.TFile(filename, 'UPDATE')
        toyTree = self.getToyTree(toysFile)
        sysTree = self.getToyTree(sysFile)
        if toyTree is not None and sysTree is not None:
            dirName = "WithToys"
            computeErrors = True
        else:
            dirName = "BeforeToys"
            computeErrors = False
        tdirectory = f.GetDirectory(dirName)
        if tdirectory==None:
            f.mkdir(dirName)
            tdirectory = f.GetDirectory(dirName)
        c = rt.TCanvas('c','c',500,400)
        rt.SetOwnership(c, False)
        rt.TH1D.SetDefaultSumw2()
        rt.TH2D.SetDefaultSumw2()
        rt.TH3D.SetDefaultSumw2()
        plotband = convertSideband(plotRegion, self.workspace,
                self.x, self.y, self.z)

        # Get the 1D and 2D histograms for each b-tag bin
        h_data_nBtagRsqMR, h_nBtagRsqMR = self.get3DFitHistos(plotband)
        h_data_MR,h_data_Rsq,h_data_RsqMR = make3DHistProjections(
                h_data_nBtagRsqMR)
        h_MR,h_Rsq,h_RsqMR = make3DHistProjections(h_nBtagRsqMR)
        if computeErrors:
            h_MR = getErrors1D(h_MR,h_data_MR,sysTree,options,"x",0,
                    len(self.x)-1,0,len(self.y)-1,0,len(self.z)-1,
                    self.x,self.y,self.z)
            h_Rsq = getErrors1D(h_Rsq,h_data_Rsq,sysTree,options,"y",0,
                    len(self.x)-1,0,len(self.y)-1,0,len(self.z)-1,
                    self.x,self.y,self.z)
            h_RsqMR = getErrors2D(h_RsqMR,h_data_RsqMR,sysTree,options,"yx",0,
                    len(self.x)-1,0,len(self.y)-1,0,len(self.z)-1,
                    self.x,self.y,self.z)
            h_nBtagRsqMR = getErrors3D(h_nBtagRsqMR,h_data_nBtagRsqMR,sysTree,
                    options,"zyx",0,len(self.x)-1,0,len(self.y)-1,0,
                    len(self.z)-1,self.x,self.y,self.z)
        for h in [h_data_MR,h_data_Rsq,h_data_RsqMR,h_data_nBtagRsqMR,
                h_MR,h_Rsq,h_RsqMR,h_nBtagRsqMR]:
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
                if computeErrors:
                    h_RsqMR_components[-1] = getErrors2D(h_RsqMR_components[-1],
                            h_data_RsqMR_components[-1],sysTree,options,
                            "yx",0,len(self.x)-1,0,len(self.y)-1,k,k,
                            self.x,self.y,self.z)
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

        # Create the nsigma histograms
        h_RsqMR_statnsigma_components = []
        h_RsqMR_nsigma_components = []
        if len(self.z)>1:
            for k in range(1,len(self.z)):
                h_RsqMR_statnsigma_components.append( 
                        getStatNSigmaHist(h_RsqMR_components[k-1], 
                            h_data_RsqMR_components[k-1],
                            "h_RsqMR_statnsigma_%ibtag"%self.z[k-1]))
                h_RsqMR_nsigma_btag = h_RsqMR.Clone(
                        "h_RsqMR_nsigma_%ibtag"%self.z[k-1])
                if computeErrors:
                    h_RsqMR_nsigma_btag = getNsigma2D(h_RsqMR_nsigma_btag,
                            h_data_RsqMR_components[k-1],toyTree,options,
                            "yx",0,len(self.x)-1,0,len(self.y)-1,k,k,
                            self.x,self.y,self.z)
                    self.blindSideband(h_RsqMR_nsigma_btag, plotRegion)
                    h_RsqMR_nsigma_components.append(h_RsqMR_nsigma_btag)

        # Print everything out to pdf files
        btagLabel = getBtagLabel(self.z)
        lumiLabel = "%.1f fb^{-1} (13 TeV)" % (self.analysis.lumi/1000.)
        boxLabel = "razor %s %s %s Fit" % (self.analysis.region,
                btagLabel,self.fitRegion.replace('LowMR,LowRsq','Sideband'))
        plotLabel = ""
        eventsLabel = "Events"
        sidebandFit = self.getSidebandMax()
        for h in h_RsqMR_components:
            tdirectory.cd()
            h.Write()
        print1DProj(c,tdirectory,h_MR,h_data_MR,
                self.dirname+"/h_MR_%s.pdf"%self.analysis.region,"M_{R} [GeV]",
                eventsLabel,lumiLabel,boxLabel,plotLabel,self.isData,False,None,
                None,h_MR_components,h_colors,h_labels)
        print1DProj(c,tdirectory,h_Rsq,h_data_Rsq,
                self.dirname+"/h_Rsq_%s.pdf"%self.analysis.region,"R^{2}",
                eventsLabel,lumiLabel,boxLabel,plotLabel,self.isData,False,None,
                None,h_Rsq_components,h_colors,h_labels)
        if len(self.z)>2:
            for k in range(0,len(self.z)-1):
                newBoxLabel = "razor %s %s %s Fit"%(self.analysis.region,
                    h_labels[k],self.fitRegion.replace('LowMR,LowRsq','Sideband'))
                if computeErrors:
                    print1DProjNs(c,tdirectory,h_th1x_components[k],
                            h_data_th1x_components[k],
                            h_RsqMR_nsigma_components[k],
                            self.dirname+"/h_th1x_ns_%ibtag_%s.pdf"%(
                                self.z[k],self.analysis.region),
                            "Bin Number",eventsLabel,lumiLabel,newBoxLabel,
                            plotLabel,self.isData,False,options,cfg=self.config) 
                    print2DResiduals(c,tdirectory,h_RsqMR_nsigma_components[k],
                            self.dirname+"/h_RsqMR_nsigma_log_%ibtag_%s.pdf"%(
                            self.z[k],self.analysis.region),"M_{R} [GeV]", 
                            "R^{2}","Stat.+Sys. n#sigma",lumiLabel,newBoxLabel,
                            plotLabel, self.x,self.y,self.isData,sidebandFit,
                            False,options)
                else:
                    print1DProj(c,tdirectory,h_th1x_components[k],
                        h_data_th1x_components[k], 
                        self.dirname+"/h_th1x_%ibtag_%s.pdf"%(self.z[k],
                            self.analysis.region),"Bin Number",
                            eventsLabel,lumiLabel,newBoxLabel,plotLabel,self.isData,
                            False, options)
                    print2DResiduals(c,tdirectory,h_RsqMR_statnsigma_components[k],
                        self.dirname+"/h_RsqMR_statnsigma_log_%ibtag_%s.pdf"%(self.z[k],
                        self.analysis.region), "M_{R} [GeV]", "R^{2}", 
                        "Stat. n#sigma (Data - Fit)/sqrt(Fit)",
                        lumiLabel,newBoxLabel,plotLabel,self.x,self.y,self.isData,
                        sidebandFit,False,options)
        f.Close()

    def doFitSequence(self, load=False, doFit=True, plot=True, unblind=False,
            runToys=False, plotToys=False):
        """Performs all steps needed to build and fit the dataset"""
        if load:
            self.loadWorkspace()
        else:
            self.initDataset()
            self.initBinnedDataset()
        if doFit:
            self.fit()
            self.plotCorrelationMatrix()
        if runToys:
            self.runToys(sysPlusStat=True)
            self.runToys(sysPlusStat=False)
        if plot:
            toysFile = None
            sysFile = None
            if runToys or plotToys:
                toysFile = self.toysFile
                sysFile = self.sysFile
            self.plot(unblind=unblind, toysFile=toysFile,
                    sysFile=sysFile)
        self.writeWorkspace()
