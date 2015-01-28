from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config


def initializeWorkspace(w,cfg):
    variables = cfg.getVariablesRange(box,"variables",w)
    parameters = cfg.getVariables(box, "parameters")
    paramNames = []
    for parameter in parameters:
        w.factory(parameter)
        paramName = parameter.split('[')[0]
        if paramName.find("Cut")==-1 and paramName.find("Ntot")==-1:
            paramNames.append(paramName)
            w.var(paramName).setConstant(False)
        else:
            if paramName.find("Ntot")==-1:
                w.var(paramName).setConstant(True)
            else:
                w.var(paramName).setConstant(False)
    pdfs = cfg.getPdfs(box,"pdfs",w)
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store fit results")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    box = options.box


    for f in args:
        if f.lower().endswith('.root'):
            rootFile = rt.TFile(f)
            workspace = rootFile.Get("w"+box)
            data = workspace.data('RMRTree')
            
    w = rt.RooWorkspace("w"+box)
    initializeWorkspace(w,cfg)
    
    rootTools.Utils.importToWS(w,data)
    
    w.Print('v')

    pdf = w.pdf('extRazorPdf')

    fitResult = pdf.fitTo(data,rt.RooFit.Save())
    fitResult.Print('v')
    
    mr = w.var('MR')
    rsq = w.var('Rsq')

    c = rt.TCanvas("c","c",600,400)
    c.SetLogy()
    mrFrame = mr.frame(400,2000,50)
    mrFrame.SetTitle("")
    mrFrame.SetXTitle("M_{R}")
    rsqFrame = rsq.frame(0.25,1.2,50)
    rsqFrame.SetTitle("")
    rsqFrame.SetXTitle("R^{2}")
    
    rootTools.Utils.importToWS(w,fitResult)
    
    def plot1d(data,pdf,var,frame,c):
        data.plotOn(frame,rt.RooFit.Name("Data"),rt.RooFit.Invisible())
        pdf.plotOn(frame,rt.RooFit.VisualizeError(fitResult,0.4),rt.RooFit.FillColor(rt.kBlue-10))
        pdf.plotOn(frame,rt.RooFit.Name("Total"),rt.RooFit.FillColor(rt.kBlue-10))
        pdf.plotOn(frame,rt.RooFit.Name("TTj1b"),rt.RooFit.Components('razor3dPdf_TTj1b'),rt.RooFit.LineColor(rt.kViolet),rt.RooFit.LineStyle(rt.kDashed))
        pdf.plotOn(frame,rt.RooFit.Name("TTj2b"),rt.RooFit.Components('razor3dPdf_TTj2b'),rt.RooFit.LineColor(rt.kRed),rt.RooFit.LineStyle(rt.kDashed))
        pdf.plotOn(frame,rt.RooFit.Name("TTj3b"),rt.RooFit.Components('razor3dPdf_TTj3b'),rt.RooFit.LineColor(rt.kGreen),rt.RooFit.LineStyle(rt.kDashed))
        data.plotOn(frame,rt.RooFit.Name("Data"))
        frame.SetMaximum(1000)
        frame.SetMinimum(0.5)
        frame.Draw()

        l = rt.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.05)
        l.SetTextFont(42)
        l.SetNDC()
        l.DrawLatex(0.15,0.85,"CMS Simulation 4 fb^{-1} (13 TeV)")
        l.DrawLatex(0.15,0.80,"Razor %s Box"%box)
        leg = rt.TLegend(0.7,0.59,0.89,0.88)
        leg.SetTextFont(42)
        leg.SetFillColor(rt.kWhite)
        leg.SetLineColor(rt.kWhite)
        leg.AddEntry(frame.findObject("Data"),"Sim Data","pe")
        leg.AddEntry(frame.findObject("Total"),"Total","lf")
        leg.AddEntry(frame.findObject("TTj1b"),"1b-tag","l")
        leg.AddEntry(frame.findObject("TTj2b"),"2b-tag","l")
        leg.AddEntry(frame.findObject("TTj3b"),"3b-tag","l")
        leg.Draw()
    
        c.Print(options.outDir+"/RooPlot_"+var.GetName()+"_"+box+".pdf")
        c.Print(options.outDir+"/RooPlot_"+var.GetName()+"_"+box+".C")

    plot1d(data,pdf,mr,mrFrame,c)
    plot1d(data,pdf,rsq,rsqFrame,c)

    inFiles = [f for f in args if f.lower().endswith('.root')]
            
    if len(inFiles)==1:
        outFile = inFiles[0].split('/')[-1].replace('RazorAnalysis','FitResult')

    outFile = rt.TFile.Open(options.outDir+"/"+outFile,'recreate')
    outFile.cd()
    w.Write()
    outFile.Close()
