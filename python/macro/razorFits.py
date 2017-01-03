import ROOT as rt 

# Local imports
from razorAnalysis import Analysis
import macro
from framework import Config
from DustinTuple2RooDataSet import convertTree2Dataset, initializeWorkspace

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
    """

    def __init__(self, box, tag="Razor2016_MoriondRereco", 
            weights=None, isData=True, configFile="config/run2_2016.config"):
        self.weights = weights
        self.isData = isData
        self.analysis = Analysis(box, tag=tag)
        self.filename = self.getFilename()
        self.config = Config.Config(configFile)
        self.workspace = rt.RooWorkspace("w"+box)
        initializeWorkspace(self.workspace, self.config, box)

    def getFilename(self):
        """Builds the name of the file to use for input/output"""
        strs = ["RazorFitInstance",self.analysis.tag,self.analysis.region]
        if not self.isData:
            strs.append("SMCocktail")
            if self.weights is not None:
                for proc,weight in self.weights.iteritems():
                    strs.append(("%s_%.2f"%(proc,weight)).replace(".","p"))
        return "_".join(strs)+".root"

    def initDataset(self):
        """Loads the ntuple files and creates the dataset for fitting"""
        if self.isData:
            filenames = { "Data":self.analysis.filenames["Data"] }
        else:
            filenames = { proc:self.analysis.filenames[proc] for proc in 
                    self.analysis.filenames if proc != "Data" and proc != "QCD" }
        files = { proc:rt.TFile.Open(filenames[proc]) for proc in filenames }
        trees = macro.makeTreeDict(files, "RazorInclusive")
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
        # Must use getattr to access the import method of the underlying 
        # ROOT object, and provide a null RooCmdArg to disambiguate
        # between the different overloaded import() methods.
        getattr(self.workspace,'import')(combinedData, rt.RooCmdArg())

    def loadWorkspace(self, filename=None):
        """Loads the workspace from file.  If no filename is given,
            the name returned by get_filename() will be used"""
        if filename is None:
            filename = self.filename
        f = rt.TFile.Open(filename)
        wName = "w"+self.box
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
        """Calls the fitting sequence"""
        raise NotImplementedError
