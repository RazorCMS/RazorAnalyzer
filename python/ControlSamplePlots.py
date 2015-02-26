#Plotting script to compare kinematic distributions in Z->nunu with those obtained from the DY+jets, W+jets, and photon+jets control samples
import sys, os
import string
import ROOT as rt

def compareControlDistributions(trees, quantitiesToPlot, titles, xmins, xmaxs, nBins, effHistogram = None, accHistogram = None, reweighBy = None, reweighByHistogram=None):
    """Draws histograms of the desired quantity for the Z->invisible samples, the DY samples, the W samples, and the Gamma samples, applies formatting, and prints the plot as a PDF.  Supply reweighByHistogram (it should be a dict of histograms returned by a different invocation of this function) to reweigh the measured distributions according to the given histograms, using the ZJets distribution as reference."""
    if reweighBy is not None and reweighByHistogram is None:
        print("Error in compareControlDistributions: reweighing expression provided without accompanying histograms to reweigh by!")
        return
    elif reweighBy is None and reweighByHistogram is not None:
        print("Error in compareControlDistributions: reweighing histograms provided without accompanying variable name!")
        return

    #names of samples to draw
    sampleNames = ["ZJets", "DYJets", "WJets", "GJets"]
    branchingRatios = {"DYJets":0.06729, "WJets":0.2132, "ZJets":0.2, "GJets":1.0}

    #legend entries
    sampleTitles = {}
    sampleTitles["ZJets"] = "Z -> #nu#nu + Jets"
    sampleTitles["DYJets"] = "DY -> ll + Jets"
    sampleTitles["WJets"] = "W -> l#nu + Jets"
    sampleTitles["GJets"] = "#gamma + Jets"

    #postfixes for each sample
    postfixes = {}
    postfixes["ZJets"] = ""
    postfixes["DYJets"] = "_noZ"
    postfixes["WJets"] = "_noW"
    postfixes["GJets"] = "_noPho"

    #efficiency expressions 
    effPtExps = {}
    effEtaExps = {}
    if effHistogram is not None:
        effPtExps["DYJets"] = "(nLooseMuons == 2)*(leadingMuonPt) + (nLooseElectrons == 2)*(leadingElectronPt)"
        effEtaExps["DYJets"] = "(nLooseMuons == 2)*abs(leadingMuonEta) + (nLooseElectrons == 2)*abs(leadingElectronEta)"
        effPtExps["WJets"] = "(nTightMuons == 1)*(leadingMuonPt) + (nTightElectrons == 1)*(leadingElectronPt)"
        effEtaExps["WJets"] = "(nTightMuons == 1)*abs(leadingMuonEta) + (nTightElectrons == 1)*abs(leadingElectronEta)"
        effPtExps["GJets"] = "(nSelectedPhotons == 1)*(leadingPhotonPt)"
        effEtaExps["GJets"] = "(nSelectedPhotons == 1)*abs(leadingPhotonEta)"
        effPtExps2 = effPtExps["DYJets"].replace("leading", "subleading") #for DY+Jets
        effEtaExps2 = effEtaExps["DYJets"].replace("leading", "subleading") #for DY+Jets

    #acceptance expressions (pt and eta of the Z boson, used to correct for lepton acceptance in DY+Jets
    accPtExps = {}
    accEtaExps = {}
    if accHistogram is not None:
        accPtExps["DYJets"] = "recoZpt"
        accEtaExps["DYJets"] = "recoZeta"

    #expression to reweigh by 
    reweighExps = {}
    if reweighBy is not None:
        for s in sampleNames: reweighExps[s] = reweighBy+postfixes[s]

    #baseline selection for each sample
    conditions = {}
    conditions["ZJets"] = "nLooseMuons == 0 && nLooseElectrons == 0 && MR > 300 && Rsq > 0.15"
    conditions["DYJets"] = "MR_noZ > 300 && Rsq_noZ > 0.15"
    conditions["WJets"] = "MR_noW > 300 && Rsq_noW > 0.15"
    conditions["GJets"] = "MR_noPho > 300 && Rsq_noPho > 0.15 && nSelectedPhotons >= 1"

    #histogram colors
    colors = {}
    colors["ZJets"] = rt.kCyan+3
    colors["DYJets"] = rt.kAzure
    colors["WJets"] = rt.kOrange+10
    colors["GJets"] = rt.kTeal+10

    #histograms, canvas and legend
    histos = {}
    for name in sampleNames: 
        histos[name] = {}
        for i, q in enumerate(quantitiesToPlot): histos[name][q] = rt.TH1F(name+q+((reweighBy is not None)*("Reweigh"))+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr")), titles[i]+"; "+titles[i], nBins[i], xmins[i], xmaxs[i])
    c = rt.TCanvas("c", "c", 800, 600)
    leg = rt.TLegend(0.7, 0.7, 0.9, 0.9)
    legNoGamma = rt.TLegend(0.7, 0.7, 0.9, 0.9)
    legNoGammaNoW = rt.TLegend(0.7, 0.7, 0.9, 0.9)

    #loop over samples and fill histograms
    for treeName in trees:
        for name in sampleNames:
            if name in treeName: #fill appropriate histogram, with weight
                if reweighByHistogram is None and effHistogram is None and accHistogram is None: #reweigh only by sample cross section
                    print("Filling "+name+" from tree "+treeName+" (weighing by sample cross section)")
                    for q in quantitiesToPlot:
                        trees[treeName].Draw(q+postfixes[name]+">>+"+histos[name][q].GetName(), "weight*("+conditions[name]+")*1.0/"+str(branchingRatios[name]))
                else: #reweigh by sample cross section and according to reweighByHistogram and/or effHistogram
                    print("Filling "+name+" from tree "+treeName+" (weighing by sample cross section "+((reweighBy is not None)*("and "+str(reweighBy)))+((effHistogram is not None)*(" and selection efficiency"))+")")
                    form = {}
                    for q in quantitiesToPlot:
                        form[q] = rt.TTreeFormula(name+q+"Formula", q+postfixes[name], trees[treeName])
                        form[q].GetNdata()
                    reweighForm = None
                    if reweighByHistogram is not None: 
                        reweighForm = rt.TTreeFormula(name+"ReweighFormula", reweighExps[name], trees[treeName])
                        reweighForm.GetNdata()
                    #formulas for reweighing based on pt and eta
                    effPtForm = None
                    effEtaForm = None
                    effPtForm2 = None
                    effEtaForm2 = None
                    accPtForm = None
                    accEtaForm = None
                    eventHasElectrons = None
                    if effHistogram is not None and name != "ZJets":
                        effPtForm = rt.TTreeFormula(name+"PtEfficiencyFormula", effPtExps[name], trees[treeName])
                        effEtaForm = rt.TTreeFormula(name+"EtaEfficiencyFormula", effEtaExps[name], trees[treeName])
                        effPtForm.GetNdata()
                        effEtaForm.GetNdata()
                        if name == "DYJets": 
                            effPtForm2 = rt.TTreeFormula(name+"PtEfficiencyFormula2", effPtExps2, trees[treeName]) #this one corrects for the efficiency of the second lepton in DY+Jets events
                            effEtaForm2 = rt.TTreeFormula(name+"EtaEfficiencyFormula2", effEtaExps2, trees[treeName]) #this one corrects for the efficiency of the second lepton in DY+Jets events
                            effPtForm2.GetNdata()
                            effEtaForm2.GetNdata()
                    if accHistogram is not None and name == "DYJets":
                        accPtForm = rt.TTreeFormula(name+"ZPtFormula", accPtExps[name], trees[treeName])
                        accEtaForm = rt.TTreeFormula(name+"ZEtaFormula", accEtaExps[name], trees[treeName])
                        accPtForm.GetNdata()
                        accEtaForm.GetNdata()

                    #determine if each event has electrons or muons
                    if name == "WJets":
                        eventHasElectrons = rt.TTreeFormula(name+"HasElectrons", "(nTightElectrons == 1)", trees[treeName]) #false for muons, true for electrons
                        eventHasElectrons.GetNdata()
                    elif name == "DYJets":
                        eventHasElectrons = rt.TTreeFormula(name+"HasElectrons", "(nLooseElectrons == 2)", trees[treeName]) #false for muons, true for electrons
                        eventHasElectrons.GetNdata()

                    conditionsForm = rt.TTreeFormula(name+"ConditionsFormula", conditions[name], trees[treeName])
                    conditionsForm.GetNdata()
                    nEvents = trees[treeName].GetEntries()
                    for n, event in enumerate(trees[treeName]):
                        #count events
                        if n % 1000000 == 0: print("Processing event "+str(n)+" of "+str(nEvents))
                        #check if the event passes
                        passesCondition = conditionsForm.EvalInstance()
                        if not passesCondition: continue
                        #value of the quantity we are filling
                        exprValue = {}
                        for q in quantitiesToPlot: exprValue[q] = form[q].EvalInstance()
                        #initialize weight
                        weight = event.weight
                        weight /= branchingRatios[name] #correct for branching fraction of Z to ele/mu, W to ele/mu, Z to nu nu
                        #value of the quantity whose distribution should be made to match that of Z+jets
                        if reweighForm is not None: #compute contribution to the weight from the reweighing histogram
                            reweighExprValue = reweighForm.EvalInstance() 
                            reweighRatio = reweighByHistogram["ZJets"].GetBinContent(reweighByHistogram["ZJets"].FindBin(reweighExprValue))*1.0/reweighByHistogram[name].GetBinContent(reweighByHistogram[name].FindBin(reweighExprValue))
                            weight = weight*reweighRatio
                        #get appropriate pt's and etas for reweighing
                        if effPtForm is not None: 
                            effPtExprValue = effPtForm.EvalInstance()
                        if effEtaForm is not None:
                            effEtaExprValue = effEtaForm.EvalInstance()
                        if effPtForm2 is not None: 
                            effPtExprValue2 = effPtForm2.EvalInstance()
                        if effEtaForm2 is not None: 
                            effEtaExprValue2 = effEtaForm2.EvalInstance()
                        if accPtForm is not None:
                            accPtExprValue = accPtForm.EvalInstance()
                        if accEtaForm is not None:
                            accEtaExprValue = accEtaForm.EvalInstance()

                        #reweigh according to selection efficiency
                        if effPtForm is not None and effEtaForm is not None and name != "ZJets":
                            #get the correct histogram for the efficiency
                            if name == "GJets":
                                thisEffHisto = effHistogram["Photon"]
                            elif eventHasElectrons and name == "WJets":
                                thisEffHisto = effHistogram["ElectronTight"]
                            elif eventHasElectrons.EvalInstance(): #DYJets
                                thisEffHisto = effHistogram["Electron"]
                            elif name == "WJets":
                                thisEffHisto = effHistogram["MuonTight"]
                            else: #DYJets
                                thisEffHisto = effHistogram["Muon"]

                            ptMax = 1000
                            if effPtExprValue > ptMax: effPtExprValue = ptMax - 1 #if pt is higher than the maximum of the histogram, use the last pt bin
                            if effPtExprValue2 > ptMax: effPtExprValue2 = ptMax - 1
                            efficiency = thisEffHisto.GetBinContent(thisEffHisto.FindBin(effPtExprValue, effEtaExprValue))
                            if name == "DYJets": efficiency *= thisEffHisto.GetBinContent(thisEffHisto.FindBin(effPtExprValue2, effEtaExprValue2))

                            if efficiency > 0.001:
                                weight = weight/efficiency
                            else: 
                                weight = 0 #unwanted event -- looking into why these happen

                        #reweigh according to acceptance
                        if accPtForm is not None and accEtaForm is not None and name == "DYJets":
                            if eventHasElectrons.EvalInstance():
                                thisAccHisto = accHistogram["Electron"]
                            else:
                                thisAccHisto = accHistogram["Muon"]

                            ptMax = 1000
                            etaMax = 6.0
                            if accPtExprValue > ptMax: effPtExprValue = ptMax - 1 
                            if abs(accEtaExprValue) > etaMax: 
                                print("Warning: Z boson pt is outside of histogram range: "+str(accEtaExprValue))
                                accEtaExprValue = etaMax - 0.1
                            accEfficiency = thisAccHisto.GetBinContent(thisAccHisto.FindBin(accPtExprValue, accEtaExprValue))
                            if accEfficiency > 0.001:
                                weight = weight/accEfficiency
                            else:
                                weight = 0

                        #multiply weight by (value of reweighing histogram in Z+jets)/(value of reweighing histogram in the given sample)
                        for q in quantitiesToPlot: histos[name][q].Fill(exprValue[q], weight)

    #format histos
    for name in sampleNames:
        for q in quantitiesToPlot:
            histos[name][q].SetLineWidth(3)
            histos[name][q].SetStats(0)
            histos[name][q].SetLineColor(colors[name])

    #add histos to legend
    for name in sampleNames: 
        #add to the legend (the itervalues.next expression just gives you one histogram from the dict
        leg.AddEntry(histos[name].itervalues().next(), sampleTitles[name])
        if name != "GJets": legNoGamma.AddEntry(histos[name].itervalues().next(), sampleTitles[name])
        if name != "GJets" and name != "WJets": legNoGammaNoW.AddEntry(histos[name].itervalues().next(), sampleTitles[name])

    #fill stacked histogram; fill a second stacked histogram with normalized versions of each histogram
    #fill a third stacked histogram with each histogram normalized to Z+Jets
    for i, q in enumerate(quantitiesToPlot):
        stack = rt.THStack("stack", titles[i]) #for drawing histograms
        normstack = rt.THStack("normstack", titles[i]) #for drawing normalized histograms
        ratiostack = rt.THStack("ratiostack", titles[i]+" ratio with Z->#nu#nu")
        ratiostackNoGamma = rt.THStack("ratiostackNoGamma", titles[i]+" ratio with Z->#nu#nu")
        ratiostackNoGammaNoW = rt.THStack("ratiostackNoGammaNoW", titles[i]+" ratio with Z->#nu#nu")
        normHistos = {}
        ratioHistos = {}
        for name in sampleNames:
            stack.Add(histos[name][q])
            normHistos[name] = histos[name][q].Clone()
            normHistos[name].Scale(1.0/normHistos[name].Integral())
            normstack.Add(normHistos[name])
            ratioHistos[name] = histos[name][q].Clone()
            ratioHistos[name].Divide(histos["ZJets"][q])
            ratioHistos[name].SetTitle("Ratio with Z->#nu#nu")
            ratiostack.Add(ratioHistos[name])
            if name != "GJets": ratiostackNoGamma.Add(ratioHistos[name])
            if name != "GJets" and name != "WJets": ratiostackNoGammaNoW.Add(ratioHistos[name])

        #draw and print all histograms
        stack.Draw("nostack,elp")
        stack.GetXaxis().SetTitle(str(titles[i]))
        leg.Draw()
        for logScale in [True, False]:
            c.SetLogy(logScale)
            c.Print("controlSample"+q+(logScale*"Log")+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+".pdf")
            c.Print("controlSample"+q+(logScale*"Log")+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+".root")

        normstack.Draw("nostack,elp")
        normstack.GetXaxis().SetTitle(str(titles[i]))
        leg.Draw()
        for logScale in [True, False]:
            c.SetLogy(logScale)
            c.Print("controlSample"+q+"Normalized"+(logScale*"Log")+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+".pdf")
            c.Print("controlSample"+q+"Normalized"+(logScale*"Log")+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+".root")

        ratiostack.Draw("nostack,elp")
        ratiostack.GetXaxis().SetTitle(str(titles[i]))
        leg.Draw()
        c.SetLogy(False)
        c.Print("controlSample"+q+"Ratio"+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+".pdf")
        c.Print("controlSample"+q+"Ratio"+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+".root")

        ratiostackNoGamma.Draw("nostack,elp")
        ratiostackNoGamma.GetXaxis().SetTitle(str(titles[i]))
        legNoGamma.Draw()
        c.SetLogy(False)
        c.Print("controlSample"+q+"Ratio"+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+"NoGamma.pdf")
        c.Print("controlSample"+q+"Ratio"+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+"NoGamma.root")

        ratiostackNoGammaNoW.Draw("nostack,elp")
        ratiostackNoGammaNoW.GetXaxis().SetTitle(str(titles[i]))
        legNoGammaNoW.Draw()
        c.SetLogy(False)
        c.Print("controlSample"+q+"Ratio"+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+"NoGammaNoW.pdf")
        c.Print("controlSample"+q+"Ratio"+((effHistogram is not None)*("EffCorr"))+((accHistogram is not None)*("AccCorr"))+((reweighBy is not None)*("ReweighBy"+str(reweighBy)))+"NoGammaNoW.root")

    #return the histograms
    return histos

##begin main program

#set ROOT to batch mode
rt.gROOT.SetBatch()

#load the efficiency histogram file; load the acceptance and efficiency histograms

efficiencyFile = rt.TFile("Phys14LeptonPhotonEfficiencyNoteIncorrectErrors.root")
effHistos = {}
accHistos = {}
for pType in ["Muon", "Electron", "Photon"]:
    effHistos[pType] = efficiencyFile.Get(pType+"Efficiency")
    if pType != "Photon": accHistos[pType] = efficiencyFile.Get(pType+"Acceptance")
effHistos["MuonTight"] = efficiencyFile.Get("MuonEfficiencyTight")
effHistos["ElectronTight"] = efficiencyFile.Get("ElectronEfficiencyTight")

#load TFiles and initialize trees
datanames = [
             "DYJets100", "DYJets200", "DYJets400", "DYJets600", 
             "ZJets100", "ZJets200", "ZJets400", "ZJets600", 
             "WJets100", "WJets200", "WJets400", "WJets600", 
             "GJets100", "GJets200", "GJets400", "GJets600"
             ]
prefix = "control"
postfix = "_4000pb_weighted.root"
files = {}
for sample in datanames: files[sample] = rt.TFile(prefix+sample+postfix)
trees = {}
treeName = "RazorInclusive"
for sample in datanames: trees[sample] = files[sample].Get(treeName)

#file for saving histograms
outfile = rt.TFile("controlSampleHistograms.root", "recreate")
outfile.cd()

#plots: MET, Rsq, MR
results = {}
#results["metHistos"] = compareControlDistributions(trees, ["met"], ["MET (GeV)"], [0.0], [1000], [40])
reweighingHistos = {}
#for key in results["metHistos"]:
#    reweighingHistos[key] = results["metHistos"][key]["met"]

#without reweighing
results["histosNoReweighing"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 4000, 15], [25, 25, 25, 25, 15])
#reweigh by lepton/photon selection efficiency
results["histosEffCorr"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 4000, 15], [25, 25, 25, 25, 15], effHistogram = effHistos)
#correct for lepton acceptance
results["histosEffCorrAccCorr"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 4000, 15], [25, 25, 25, 25, 15], effHistogram = effHistos, accHistogram = accHistos)
#after reweighing by MET
for name in results["histosEffCorrAccCorr"]: reweighingHistos[name] = results["histosEffCorrAccCorr"][name]["met"] #set up reweighing histogram using the latest MET distribution
results["histosEffCorrAccCorrMetCorr"] = compareControlDistributions(trees, ["met", "Rsq", "MR", "HT", "numJets"], ["MET (GeV)", "R^{2}", "M_{R} (GeV)", "HT (GeV)", "Number of jets"], [0.0, 0.0, 0.0, 0.0, 0], [1000, 1, 1500, 4000, 15], [25, 25, 25, 25, 15], effHistogram = effHistos, accHistogram = accHistos, reweighBy = "met", reweighByHistogram = reweighingHistos)

#save histograms to ROOT file
for result in results:
    for key in results[result]:
        for q in results[result][key]:
            results[result][key][q].Write()
