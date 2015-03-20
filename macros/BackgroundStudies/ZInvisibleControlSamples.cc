//Macro to predict the MR and Rsq distributions of the Z->nu nu background in the razor search using Z->mu mu, W->mu nu, and Gamma+Jets events

#include <iostream>
#include <map>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

using namespace std;

void ZInvisibleControlSamples(){
    //decide to reweigh by MET or by MR and Rsq
    //(reweighByRazor = false to reweigh by MET, true to reweigh by MR and Rsq)
    bool reweighByRazor = false; 

    map<string, string> suffixes;
    suffixes["DYJets"] = "_noZ";
    suffixes["WJets"] = "_noW";
    suffixes["GJets"] = "_noGamma";
    map<string, float> normFactors; //normalize DYJets and WJets to Z->nu nu cross section
    normFactors["DYJets"] = 0.03366/0.2;
    normFactors["WJets"] = 2363.02/489.503*10.57/32.57;
    //get input files -- assumes one TFile for each process, with weights for different HT bins 
    map<string, TFile*> mcfiles;
    map<string, TFile*> datafiles;
    mcfiles["DYJets"] = new TFile("DYJets.root");
    mcfiles["WJets"] = new TFile("WJets.root");
    //mcfiles["GJets"] = new TFile("GJets.root");
    datafiles["DYJets"] = new TFile("DoubleMu.root");
    datafiles["WJets"] = new TFile("SingleMu.root");
    //datafiles["GJets"] = new TFile("Photon.root");
    //get trees and set branches
    map<string, TTree*> mctrees;
    map<string, TTree*> datatrees;
    map<string, float> mets;
    map<string, float> mrs;
    map<string, float> rsqs;
    float weight;
    float leadingMuonPt, leadingMuonEta, leadingPhotonPt, leadingPhotonEta, recoZpt, recoZeta, subleadingMuonPt, subleadingMuonEta;
    for(auto &file : mcfiles){
        mets[file.first] = 0.;
        mrs[file.first] = 0.;
        rsqs[file.first] = 0.;
        mctrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        mctrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        mctrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
        mctrees[file.first]->SetBranchAddress("weight", &weight);
        mctrees[file.first]->SetBranchAddress("leadingMuonPt", &leadingMuonPt);
        mctrees[file.first]->SetBranchAddress("leadingMuonEta", &leadingMuonEta);
        mctrees[file.first]->SetBranchAddress("leadingPhotonPt", &leadingPhotonPt);
        mctrees[file.first]->SetBranchAddress("leadingPhotonEta", &leadingPhotonEta);
        mctrees[file.first]->SetBranchAddress("recoZpt", &recoZpt);
        mctrees[file.first]->SetBranchAddress("recoZeta", &recoZeta);
    }
    for(auto &file : datafiles){
        datatrees[file.first] = (TTree*)file.second->Get("RazorInclusive");
        datatrees[file.first]->SetBranchAddress(Form("met%s", suffixes[file.first].c_str()), &mets[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("MR%s", suffixes[file.first].c_str()), &mrs[file.first]);
        datatrees[file.first]->SetBranchAddress(Form("Rsq%s", suffixes[file.first].c_str()), &rsqs[file.first]);
    }

    //load efficiency/acceptance histograms
    TFile effFile("Run1LeptonPhotonEfficiencyNoteIncorrectErrors.root");
    TH2F muonLooseEffHisto = *(TH2F *)effFile.Get("MuonEfficiency");
    TH2F muonTightEffHisto = *(TH2F *)effFile.Get("MuonEfficiencyTight");
    //TH2F photonEffHisto = *(TH2F *)effFile.Get("PhotonEfficiency");
    //TH2F zAccHisto = *(TH2F *)effFile.Get("ZAcceptance");
    //TH2F wAccHisto = *(TH2F *)effFile.Get("WAcceptance");
    
    //Step 1: Get the distributions to reweigh by: MET, MR, Rsq
    map<string, TH1F> metHistosForReweighing;
    map<string, TH2F> razorHistosForReweighing;
    for(auto &tree : mctrees){
        cout << "Filling MC histograms: " << tree.first << endl;
        metHistosForReweighing[tree.first] = TH1F(Form("metmc%s", tree.first.c_str()), "MET (GeV); MET(GeV)", 40, 0., 1000);
        razorHistosForReweighing[tree.first] = TH2F(Form("razormc%s", tree.first.c_str()), "; MR (GeV); Rsq", 40, 300., 4000, 40, 0.15, 1.5);
        uint nEntries = tree.second->GetEntries();
        for(uint i = 0; i < nEntries; i++){
            tree.second->GetEntry(i); 
            float eventWeight = weight/normFactors[tree.first];
            //reweigh according to selection efficiency and acceptance
            /*if(tree.first == "GJets"){
                eventWeight *= photonEffHisto.GetBinContent(photonEffHisto.FindBin(leadingPhotonPt, leadingPhotonEta));
            }
            else if(tree.first == "WJets"){
                eventWeight *= muonTightEffHisto.GetBinContent(muonTightEffHisto.FindBin(leadingMuonPt, leadingMuonEta));
            }
            else if(tree.first == "DYJets"){
                eventWeight *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(leadingMuonPt, leadingMuonEta));
                eventWeight *= muonLooseEffHisto.GetBinContent(muonLooseEffHisto.FindBin(subleadingMuonPt, subleadingMuonEta));
                eventWeight *= zAccHisto.GetBinContent(zAccHisto.FindBin(recoZpt, recoZeta));
            }
            else{
                cout << "Error in efficiency reweighing; check the code" << endl;
            }*/
            //fill each quantity
            metHistosForReweighing[tree.first].Fill(mets[tree.first], eventWeight);
            razorHistosForReweighing[tree.first].Fill(mrs[tree.first], rsqs[tree.first], eventWeight);
        }
    }

    //Step 3: Apply the reweighing factors to data
    map<string, TH2F> razorHistosData;
    for(auto &tree : datatrees){
        cout << "Filling data histograms: " << tree.first << endl;
        razorHistosData[tree.first] = TH2F(Form("razordata%s", tree.first.c_str()), "; MR (GeV); Rsq", 25, 300., 4000, 25, 0.15, 1.5);
        uint nEntries = tree.second->GetEntries();
        for(uint i = 0; i < nEntries; i++){
            tree.second->GetEntry(i);
            float reweighFactor = 0.0;
            if(reweighByRazor){ //reweigh by MR and Rsq
                //get the factor to reweigh by
                float numerator = razorHistosForReweighing[tree.first].GetBinContent(razorHistosForReweighing[tree.first].FindBin(mrs[tree.first], rsqs[tree.first]));
                float denominator = razorHistosForReweighing["DYJets"].GetBinContent(razorHistosForReweighing["DYJets"].FindBin(mrs[tree.first], rsqs[tree.first]));
                if(denominator > 0){
                    reweighFactor = numerator / denominator;
                }
            } 
            else{ //reweigh by MET
                //get the factor to reweigh by
                float numerator = metHistosForReweighing[tree.first].GetBinContent(metHistosForReweighing[tree.first].FindBin(mets[tree.first]));
                float denominator = metHistosForReweighing["DYJets"].GetBinContent(metHistosForReweighing["DYJets"].FindBin(mets[tree.first]));
                if(denominator > 0){
                    reweighFactor = numerator / denominator;    
                }
            }
            razorHistosData[tree.first].Fill(mrs[tree.first], rsqs[tree.first], reweighFactor);
        }
        TFile outfile("controlSampleHistograms.root", "recreate");
        TCanvas c("c", "c", 800, 600);
        c.SetLogx();
        for(auto &hist : razorHistosData){
            hist.second.Draw("text");
            c.Print(Form("controlSampleHistogram%s.pdf", hist.first.c_str()));
            c.Print(Form("controlSampleHistogram%s.root", hist.first.c_str()));
            hist.second.Write();
        }
    }

}

int main(){
    ZInvisibleControlSamples();
    return 0;
}
