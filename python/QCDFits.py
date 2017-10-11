import numpy as np
import ROOT as rt


def getXbins(hist):
    """Returns array of bin edges on x-axis for hist"""
    xbins = []
    for i in range(hist.GetNbinsX()):
        xbins.append( hist.GetXaxis().GetXbins()[i] )
    xbins.append(4000)
    return np.array(xbins) 

def makeQCDHist(name, xbins):
    """Makes histogram with the given name and bin edges
        on the x-axis, and b-tags (0,1,2) on the y-axis"""
    ybins = np.array([0.0, 1.0, 2.0, 3.0])
    hist = rt.TH2F(name, name, len(xbins)-1, xbins, len(ybins)-1, ybins)
    hist.GetXaxis().SetTitle("MR (GeV)")
    hist.GetYaxis().SetTitle("Number of b-tags")
    return hist


if __name__ == '__main__':
    rt.gROOT.SetBatch()

    inFile = rt.TFile.Open("/afs/cern.ch/work/j/jlawhorn/public/Razor_Moriond2017/clean/CMSSW_7_1_5/src/RazorAnalyzer/macros/BackgroundStudies/QCD/qcdTranslationFactors_final.root")
    outFile = rt.TFile.Open("data/ScaleFactors/RazorMADD2015/RazorQCDScaleFactors_Razor2016_MoriondRereco.root", "recreate")
    c = rt.TCanvas("c","c",800,600)
    for box in ['MultiJet','DiJet']:
        testHist = inFile.Get("npf_2d_{}_0b".format(box.lower()))
        xbins = getXbins(testHist)

        slopes = makeQCDHist("QCDSlopes_{}".format(box), xbins)
        inters = makeQCDHist("QCDInters_{}".format(box), xbins)
        covars = makeQCDHist("QCDCovars_{}".format(box), xbins)

        for btags in range(3):
            inHist = inFile.Get("npf_2d_{}_{}b".format(box.lower(), btags))
            
            for ibin in range(1, len(xbins)):
                print "{} box {} btags, {} < MR < {}".format(
                        box, btags, xbins[ibin-1], xbins[ibin])
                histSlice = inHist.ProjectionY("{}_{}_bin{}".format(
                    box, btags, ibin), ibin, ibin)
                if histSlice.Integral() == 0:
                    print "Slice has no data"
                else:
                    result = histSlice.Fit("pol1", "smf")
                    params = result.GetParams()
                    errors = result.GetErrors()
                    covar   = result.GetCovarianceMatrix()(0,1)
                    inters.SetBinContent(ibin, btags+1, params[0])
                    inters.SetBinError(ibin, btags+1, errors[0])
                    slopes.SetBinContent(ibin, btags+1, params[1])
                    slopes.SetBinError(ibin, btags+1, errors[1])
                    covars.SetBinContent(ibin, btags+1, covar)
                    covars.SetBinError(ibin, 0.)
        print "Slopes"
        slopes.Print("all")
        slopes.Draw("colztext")
        c.Print("QCDSlopes_{}.pdf".format(box))
        print "Intercepts"
        inters.Print("all")
        inters.Draw("colztext")
        c.Print("QCDIntercepts_{}.pdf".format(box))
        print "Covariances"
        covars.Print("all")
        covars.Draw("colztext")
        c.Print("QCDCovariances_{}.pdf".format(box))

        outFile.cd()
        slopes.Write()
        inters.Write()
        covars.Write()
