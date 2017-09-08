import ROOT as rt
import os
from array import array 
rt.gROOT.SetBatch()

#myfile = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p8_24Apr2017_Trigger/OneLeptonFull/RazorDM_OneLeptonFull_SingleLeptonSkim_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root")
myfile = rt.TFile.Open("/eos/cms/store/group/phys_susy/razor/Run2Analysis/RunTwoRazorControlRegions/2016/V3p8_18Jul2017_Trigger_SingleElectron2016/OneLeptonFull/RazorDM_OneLeptonFull_SingleLeptonSkim_Razor2016_MoriondRereco_Data_NoDuplicates_RazorSkim_GoodLumiGolden.root")

assert(myfile)
mytree = myfile.Get("ControlSampleEvent")
MR_bin = [0.,100.,150.,200.,250.,300.,350.,400.,550.,800.]
Rsq_bin = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1,1.2]
MR_eff = rt.TEfficiency("MR_eff",";M_{R}",len(MR_bin)-1,array('d',MR_bin))
HT_eff = rt.TEfficiency("HT_eff",";HT",60,0,2500)
Rsq_eff = rt.TEfficiency("Rsq_eff",";R^{2}",len(Rsq_bin)-1,array('d',Rsq_bin))
MRRsq_eff = rt.TEfficiency("MRRssq_eff",";M_{R};R^{2}",40,0,800,40,0,1.)
MR_Rsq_gt_0p5_eff = rt.TEfficiency("MR_Rsq_gt_0p03_eff",";M_{R}",len(MR_bin)-1, array('d', MR_bin))
Rsq_MR_gt_200_eff = rt.TEfficiency("Rsq_MR_gt_200_eff",";R^{2}",len(Rsq_bin)-1, array('d', Rsq_bin))

Rsq_200_lt_MR_lt_250_eff = rt.TEfficiency("Rsq_200_lt_MR_lt_250_eff",";R^{2}",len(Rsq_bin)-1, array('d', Rsq_bin))
Rsq_250_lt_MR_lt_350_eff = rt.TEfficiency("Rsq_250_lt_MR_lt_350_eff",";R^{2}",len(Rsq_bin)-1, array('d', Rsq_bin))
Rsq_350_lt_MR_lt_450_eff = rt.TEfficiency("Rsq_350_lt_MR_lt_450_eff",";R^{2}",len(Rsq_bin)-1, array('d', Rsq_bin))
Rsq_MR_gt_450_eff = rt.TEfficiency("Rsq_MR_gt_450_eff",";R^{2}",len(Rsq_bin)-1, array('d', Rsq_bin))

MR_eff.SetStatisticOption(rt.TEfficiency.kFCP)
MR_eff.SetConfidenceLevel(0.68)

Rsq_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_eff.SetConfidenceLevel(0.68)

MRRsq_eff.SetStatisticOption(rt.TEfficiency.kFCP)
MRRsq_eff.SetConfidenceLevel(0.68)

Rsq_MR_gt_200_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_MR_gt_200_eff.SetConfidenceLevel(0.68)

Rsq_200_lt_MR_lt_250_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_200_lt_MR_lt_250_eff.SetConfidenceLevel(0.68)

Rsq_250_lt_MR_lt_350_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_250_lt_MR_lt_350_eff.SetConfidenceLevel(0.68)

Rsq_350_lt_MR_lt_450_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_350_lt_MR_lt_450_eff.SetConfidenceLevel(0.68)

Rsq_MR_gt_450_eff.SetStatisticOption(rt.TEfficiency.kFCP)
Rsq_MR_gt_450_eff.SetConfidenceLevel(0.68)

MR_Rsq_gt_0p5_eff.SetStatisticOption(rt.TEfficiency.kFCP)
MR_Rsq_gt_0p5_eff.SetConfidenceLevel(0.68)

print "NEntries = %d" %(mytree.GetEntries())
for i in range(mytree.GetEntries()/10):
    mytree.GetEntry(i)
    if (i%10000==0): print("Get entry %d" %(i))
    #if ((mytree.HLTDecision[150] or mytree.HLTDecision[151] or mytree.HLTDecision[152] or mytree.HLTDecision[153] or mytree.HLTDecision[154])):
    if (mytree.HLTDecision[27] or mytree.HLTDecision[29] or mytree.HLTDecision[34] or mytree.HLTDecision[36] or mytree.HLTDecision[37] or mytree.HLTDecision[38] or mytree.HLTDecision[39] or mytree.HLTDecision[42]) and mytree.MR > 0 and mytree.Rsq > 0: # Ele
    #if (mytree.HLTDecision[4] or mytree.HLTDecision[13] or mytree.HLTDecision[18] or mytree.HLTDecision[20] or mytree.HLTDecision[24] or mytree.HLTDecision[29] or mytree.HLTDecision[29] or mytree.HLTDecision[34] or mytree.HLTDecision[36] or mytree.HLTDecision[37] or mytree.HLTDecision[38] or mytree.HLTDecision[39] or mytree.HLTDecision[42]) and mytree.NJets80 > 1: # 1L and 2L
    #if (mytree.HLTDecision[150] or mytree.HLTDecision[151] or mytree.HLTDecision[117] or mytree.HLTDecision[118] or mytree.HLTDecision[119] or mytree.HLTDecision[127] or mytree.HLTDecision[128] or mytree.HLTDecision[129] ): #PFHT125 || PFHT200 || PFJet and Dijet 40 60 80
        if mytree.MR > 200:
            passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
            Rsq_MR_gt_200_eff.Fill(passSel, mytree.Rsq)
            if mytree.MR < 250:
                passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
                Rsq_200_lt_MR_lt_250_eff.Fill(passSel, mytree.Rsq)
        
        if mytree.MR > 250 and mytree.MR < 350:
            passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
            Rsq_250_lt_MR_lt_350_eff.Fill(passSel, mytree.Rsq)

        if mytree.MR > 350 and mytree.MR < 450:
            passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
            Rsq_350_lt_MR_lt_450_eff.Fill(passSel, mytree.Rsq)

        if mytree.MR > 450:
            passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
            Rsq_MR_gt_450_eff.Fill(passSel, mytree.Rsq)


        if mytree.Rsq > 0.5:
            passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
            MR_Rsq_gt_0p5_eff.Fill(passSel, mytree.MR)
            
        passSel = mytree.HLTDecision[172] or mytree.HLTDecision[166]
        MRRsq_eff.Fill(passSel, mytree.MR, mytree.Rsq)
        MR_eff.Fill(passSel, mytree.MR)
        Rsq_eff.Fill(passSel, mytree.Rsq)
        HT_eff.Fill(passSel, mytree.HT)
        
rt.gStyle.SetOptStat(0)
c1 = rt.TCanvas("c1","",600,600)
rt.gPad.SetGridx()
rt.gPad.SetGridy()

Rsq_MR_gt_200_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_MR_gt_200.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_MR_gt_200.C")

Rsq_MR_gt_450_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_MR_gt_450.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_MR_gt_450.C")

Rsq_200_lt_MR_lt_250_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_200_lt_MR_lt_250.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_200_lt_MR_lt_250.C")

Rsq_250_lt_MR_lt_350_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_250_lt_MR_lt_350.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_250_lt_MR_lt_350.C")

Rsq_350_lt_MR_lt_450_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_350_lt_MR_lt_450.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq_350_lt_MR_lt_450.C")

MR_Rsq_gt_0p5_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MR_Rsq_gt_0p5.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MR_Rsq_gt_0p5.C")

Rsq_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_Rsq.C")

MR_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MR.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MR.C")

HT_eff.Draw("AP")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_HT.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_HT.C")

MRRsq_eff.Draw("COLZ")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MRRsq.png")
c1.SaveAs("/eos/user/q/qnguyen/www/eff_MRRsq.C")
