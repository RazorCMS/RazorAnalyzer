#include "RazorAnalyzer.h"
#include "TLorentzVector.h"

double RazorAnalyzer::deltaPhi(double phi1, double phi2) {
  double dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  return dphi;
}

double RazorAnalyzer::deltaR(double eta1, double phi1, double eta2, double phi2) {
  double dphi = deltaPhi(phi1,phi2);
  double deta = eta1 - eta2;
  return sqrt( dphi*dphi + deta*deta);
}

TLorentzVector RazorAnalyzer::makeTLorentzVector(double pt, double eta, double phi, double energy){
    double px = pt*cos(phi);
    double py = pt*sin(phi);
    double pz = pt*sinh(eta);

    TLorentzVector vec(px, py, pz, energy);
    return vec;
}

TLorentzVector RazorAnalyzer::makeTLorentzVectorPtEtaPhiM(double pt, double eta, double phi, double mass){
    TLorentzVector vec;
    vec.SetPtEtaPhiM(pt, eta, phi, mass);
    return vec;
}

vector<TLorentzVector> RazorAnalyzer::getHemispheres(vector<TLorentzVector> jets){
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;

    //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
    int nComb = pow(2, nJets);

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }

    //step 2: choose the partition that minimizes m1^2 + m2^2
    double mMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;
    for(size_t i=0; i < possibleHem1s.size(); i++){
        double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
        if(mMin < 0 || mTemp < mMin){
            mMin = mTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }

    //return the hemispheres in decreasing order of pt
    vector<TLorentzVector> hemsOut;
    if(myHem1.Pt() > myHem2.Pt()){
        hemsOut.push_back(myHem1);
        hemsOut.push_back(myHem2);
    } else {
        hemsOut.push_back(myHem2);
        hemsOut.push_back(myHem1);
    }

    return hemsOut;
}

double RazorAnalyzer::computeMR(TLorentzVector hem1, TLorentzVector hem2){
    return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

double RazorAnalyzer::computeRsq(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
    double mR = computeMR(hem1, hem2);
    double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
    double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
    double mTR = sqrt(term1 - term2);
    return (mTR / mR) * (mTR / mR);
}

//auxiliary functions for RazorInclusive and MatchedRazorInclusive analyses
bool RazorAnalyzer::passesHadronicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 0 || Rsq < 0) passes = false;
    //temporarily disable these
    //if(MR < 400 || Rsq < 0.25) passes = false;
    //if(MR < 450 && Rsq < 0.3) passes = false;
    return passes;
}

bool RazorAnalyzer::passesLeptonicRazorBaseline(double MR, double Rsq){
    bool passes = true;
    if(MR < 0 || Rsq < 0) passes = false;
    //temporarily disable these
    // if(MR < 300 || Rsq < 0.15) passes = false;
    // if(MR < 350 && Rsq < 0.2) passes = false;
    return passes;
}

//Checks if ToSubtract matches any particles in Collection, and subtracts the momentum of ToSubtract from the closest one
int RazorAnalyzer::SubtractParticleFromCollection(TLorentzVector ToSubtract, vector<TLorentzVector>& Collection, float deltaRMatch){
    //if Collection has any elements within R<deltaRMatch of the vector ToSubtract, find the closest such element
    //otherwise, return -1
    double closestDR = -1;
    int closestDRIndex = -1;
    for(uint i = 0; i < Collection.size(); i++){
        double thisDR = Collection[i].DeltaR(ToSubtract);
        if(closestDR < 0 || thisDR < closestDR){
            closestDR = thisDR;
            closestDRIndex = i;
        }
    }
    if(closestDR < 0 || closestDR > deltaRMatch){ //if we didn't look at any objects or we didn't get one within 0.4, return -1
        return -1;
    }
    
    //subtract the momentum (magnitude) of ToSubtract from that of the closest vector in Collection

    //if we're subtracting away everything, just set the vector to 0 and return the index of the changed vector
    if(ToSubtract.P() >= Collection[closestDRIndex].P()){ 
        Collection[closestDRIndex].SetPtEtaPhiE(0.0, 0.0, 0.0, 0.0);
        return closestDRIndex;
    }
    //otherwise, scale the 4-momentum by the appropriate factor
    double scalePFactor = (Collection[closestDRIndex].P() - ToSubtract.P())/Collection[closestDRIndex].P(); // = new P / old P
    Collection[closestDRIndex].SetPxPyPzE(Collection[closestDRIndex].Px()*scalePFactor, Collection[closestDRIndex].Py()*scalePFactor, Collection[closestDRIndex].Pz()*scalePFactor, Collection[closestDRIndex].E()*scalePFactor);
    return closestDRIndex;
}
