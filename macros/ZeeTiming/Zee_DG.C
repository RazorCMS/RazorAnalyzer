
#include <iostream>
#include <fstream> 
#include <sstream>
#include <map>
#include <utility>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLatex.h"
#include <math.h>
#include <time.h>
#include <algorithm>
#include <functional>
#include "TRandom3.h"
#include "TError.h"
#include <cmath>


TF1* Fitter(TH1 *hist, bool single = 0) {
  float xmin = hist->GetMean() - 3.0*hist->GetRMS();
  float xmax = hist->GetMean() + 3.0*hist->GetRMS();
  // fitting SINGLE Gaussian
  if(single == 1) {
    hist->Fit("gaus","Q","",xmin,xmax); // Q suppresses fit results
    gStyle->SetOptFit(1);
    return hist->GetFunction("gaus"); 
  }
  // fitting DOUBLE Gaussian
  TF1 *DGfit = new TF1("DGfit", "gaus(0) + gaus(3)", xmin, xmax);
  DGfit->SetParameters(hist->GetMaximum()/2.0, 0, hist->GetRMS(),
                       hist->GetMaximum()/2.0, 0, hist->GetRMS());
  DGfit->SetParLimits(1, -1.0, 1.0 );
  DGfit->SetParLimits(4, -1.0, 1.0 );
  DGfit->SetParLimits(2, .001, 3.0 );
  DGfit->SetParLimits(5, .001, 3.0 );
  DGfit->SetParNames("A_{1}", "#mu_{1}", "#sigma_{1}", "A_{2}", "#mu_{2}", "#sigma_{2}");
  gPrintViaErrorHandler = kTRUE; // Deal with MINOS print statements in error handler (for next line to be useful)
  gErrorIgnoreLevel = kWarning; //Suppress MINOS new minimum message
  hist->Fit("DGfit","Q","",xmin,xmax); // Q suppresses fit results. Was QMLES
  gErrorIgnoreLevel = kInfo;
  gStyle->SetOptFit(1);
  return hist->GetFunction("DGfit");
  // Call output with, e.g. ->GetParameter(2)
}



void scatterPlot(float x[], float y[], float ex[], float ey[], int points, TH1F *hist[], string xTitle, string yTitle, string Title, string filename) {
  // Copy input arrays, so changes won't affect the initial arrays
  float X[points];
  float Y[points];
  float eX[points];
  float eY[points];
  std::copy(x, x+points, X);
  std::copy(y, y+points, Y);
  std::copy(ex, ex+points, eX);
  std::copy(ey, ey+points, eY);

  // Correct For Low-Event Histograms:
  int  nBad = 0;
  for (int i=0; i<points; i++) {
    if (hist[i]->GetEntries() < 1E3) {
      for (int j=i; j<(points-1); j++){
        X[j-nBad] = X[j-nBad+1]; // delete bad element and shift others left
        Y[j-nBad] = Y[j-nBad+1]; // Note last element gets duplicated but also stays where it is
        eX[j-nBad]=eX[j-nBad+1]; // This is ok, since we will plot the first (points-nBad) points
        eY[j-nBad]=eY[j-nBad+1];
      }
      nBad +=1;
    }
  }
  if (nBad == points) {cout<<"Plot could not be made. Too few statistics."<<endl; return;}

  // Make the Plot:
  TCanvas *c = new TCanvas("c","c",200,10,600,400);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.15);
  c->SetTopMargin(0.12);
  c->cd();

  TGraphErrors* ge = new TGraphErrors(points-nBad, X, Y, eX, eY);
  ge->Draw("ap");

  TAxis* Xaxis = ge->GetXaxis();
  Xaxis->SetTitle( xTitle.c_str() );
  Xaxis->SetTitleSize(0.08);
  Xaxis->SetLabelSize(0.06);
  Xaxis->SetTitleOffset(0.8);

  TAxis* Yaxis = ge->GetYaxis();
  Yaxis->SetTitle( yTitle.c_str() );
  Yaxis->SetTitleSize(0.08);
  Yaxis->SetLabelSize(0.06);
  Yaxis->SetTitleOffset(0.85);

  gStyle->SetTitleFontSize(0.1);
  ge->SetTitle( Title.c_str() );
  if (points-nBad < 25) ge->SetLineWidth(3);
  else if (points-nBad < 50) ge->SetLineWidth(2);
  else ge->SetLineWidth(1);

  gErrorIgnoreLevel = kWarning; //Suppress info about created file
  c->SaveAs( ("plots/pdf/"+filename+".pdf").c_str() );
  c->SaveAs( ("plots/png/"+filename+".png").c_str() );
  c->SaveAs( ("plots/c/"+filename+".C").c_str() );
  delete c;
  gErrorIgnoreLevel = kInfo;
}



float GetDGSigma(TF1 *fit, bool single = 0) {
  if (single==1) return fit->GetParameter(2);
  float p[6];
  for (int i=0; i<6; i++) p[i] = fit->GetParameter(i);
  return (p[0] * p[2] + p[3] * p[5]) / (p[0] + p[3]);
}



float GetDGMean(TF1 *fit, bool single = 0) {
  if (single==1) return fit->GetParameter(1);
  float p[6];
  for (int i=0; i<6; i++) p[i] = fit->GetParameter(i);
  return (p[0] * p[1] + p[3] * p[4]) / (p[0] + p[3]);
}



float GetDGSigmaError(TF1 *fit, bool single = 0) {
  if (single==1) return fit->GetParError(2);
  float p[6];
  float e[6];
  for (int i=0; i<6; i++) {
    p[i] = fit->GetParameter(i);
    e[i] = fit->GetParError(i);
  }
  return sqrt(pow(p[0]*e[2], 2) + pow(p[3]*e[5], 2) + pow( (p[2]-p[5])/(p[0]+p[3]), 2) * ( pow(p[3]*e[0], 2) + pow(p[0]*e[3], 2) ))
         / (p[0] + p[3]);
}



float GetDGMeanError(TF1 *fit, bool single = 0) {
  if (single==1) return fit->GetParError(1);
  float p[6];
  float e[6];
  for (int i=0; i<6; i++) {
    p[i] = fit->GetParameter(i);
    e[i] = fit->GetParError(i);
  }
  return sqrt(pow(p[0]*e[1], 2) + pow(p[3]*e[4], 2) + pow( (p[1]-p[4])/(p[0]+p[3]), 2) * ( pow(p[3]*e[0], 2) + pow(p[0]*e[3], 2) ))
         / (p[0] + p[3]);
}



template <class T1, class T2>
void SigMeanTGraph(TFile *file, TH1F** hists, string histName, T1 cutArray[], T2 stepSize, int nSteps, string Time, string xTitle, string outFile) {
  TF1 *Fit[nSteps];
  float  Sigma[nSteps];
  float  SigEr[nSteps];
  float   Mean[nSteps];
  float MeanEr[nSteps];
  bool singleGauss = 0;

  for (int i=0; i<nSteps; i++) {
    if ( !(hists[i]->GetEntries() > 100) ) continue;
    Fit[i] = Fitter( hists[i], singleGauss );
    file->WriteTObject(hists[i], ( histName+Form("[%d]",i) ).c_str(),"WriteDelete");
    Sigma[i] = GetDGSigma( Fit[i], singleGauss );     
    SigEr[i] = GetDGSigmaError( Fit[i], singleGauss ); 
    Mean[i]  = GetDGMean( Fit[i], singleGauss );   
    MeanEr[i]= GetDGMeanError( Fit[i], singleGauss );  
  }

  float xVals[nSteps];
  float xErr[nSteps];
  for(int i=0; i<nSteps; i++) xVals[i] = cutArray[i] + stepSize/2;
  std::fill(xErr, xErr+nSteps, float(stepSize)/2);

  scatterPlot(xVals, Sigma, xErr, SigEr, nSteps, hists, xTitle.c_str(), "Time Resolution (ns)", ("#sigma of t_{1}-t_{2} "+Time).c_str() , (outFile+"_Sigma").c_str() );
  scatterPlot(xVals, Mean, xErr, MeanEr, nSteps, hists, xTitle.c_str(), "Mean", ("Mean of t_{1}-t_{2} "+Time).c_str() , (outFile+"_Mean").c_str() );
}



TF1* SaveHist(TH1 *hist, int linearfit = 0, float xmin = 1.3, float xmax = 2.3) {
  TCanvas *c = new TCanvas("c","c",200,10,600,400);
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.15);
  c->SetTopMargin(0.12);
  c->SetRightMargin(0.12);
  c->cd();
  
  hist->Draw();

  TAxis* Xaxis = hist->GetXaxis();
  Xaxis->SetTitleSize(0.08);
  Xaxis->SetLabelSize(0.06);
  Xaxis->SetTitleOffset(0.8);

  TAxis* Yaxis = hist->GetYaxis();
  Yaxis->SetTitleSize(0.08);
  Yaxis->SetLabelSize(0.06);
  Yaxis->SetTitleOffset(0.85);

  TF1 *linfit;
  if (linearfit == 1) { 
    linfit = new TF1("linfit", "pol1", xmin, xmax);
    hist->Fit("linfit","QMES", "", xmin, xmax);//linfit","QMES", 1.3, 2.3);
    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.080);
    tex->SetTextFont(42);
    tex->SetTextColor(kRed);
    float yint = linfit->GetParameter(0);
    float slope= linfit->GetParameter(1);
    if( yint > 0 ) tex->DrawLatex( 0.15, 0.8, Form("y = %.2fx + %.2f", slope, yint) );
    else tex->DrawLatex( 0.15, 0.8, Form("y = %.2fx - %.2f", slope, -yint) );
  }

  gStyle->SetTitleFontSize(0.1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gErrorIgnoreLevel = kWarning; //Suppress info about created file
  c->SaveAs( (string("plots/pdf/")+ hist->GetName() +".pdf").c_str() );
  c->SaveAs( (string("plots/png/")+ hist->GetName() +".png").c_str() );
  c->SaveAs( (string("plots/c/")  + hist->GetName() +".C").c_str() );
  delete c;
  gErrorIgnoreLevel = kInfo;
  
  return linfit;
}



void EtaPhi1D( int etaChannels, int phiChannels, TH1F** etaHists, TH1F** phiHists, string timeTitle ) {
  TF1 *etaFit[etaChannels];
  TF1 *phiFit[phiChannels];
  float etaSigma[etaChannels];
  float etaSigEr[etaChannels];
  float phiSigma[phiChannels];
  float phiSigEr[phiChannels];
  float etaMean[etaChannels];
  float etaMeanEr[etaChannels];
  float phiMean[phiChannels];
  float phiMeanEr[phiChannels];
  bool singleGauss = 0;

  for (int i=0; i<phiChannels; i++) {
    if ( !(phiHists[i]->GetEntries() != 0) ) continue;
    phiFit[i] = Fitter( phiHists[i], singleGauss );
    phiSigma[i] = GetDGSigma( phiFit[i], singleGauss );      
    phiSigEr[i] = GetDGSigmaError( phiFit[i], singleGauss ); 
    phiMean[i]  = GetDGMean( phiFit[i], singleGauss );       
    phiMeanEr[i]= GetDGMeanError( phiFit[i], singleGauss );  
  }
  for (int i=0; i<etaChannels; i++) {
    if ( !(etaHists[i]->GetEntries() != 0) ) continue;
    etaFit[i] = Fitter( etaHists[i], singleGauss );
    etaSigma[i] = GetDGSigma( etaFit[i], singleGauss );
    etaSigEr[i] = GetDGSigmaError( etaFit[i], singleGauss );
    etaMean[i]  = GetDGMean( etaFit[i], singleGauss );
    etaMeanEr[i]= GetDGMeanError( etaFit[i], singleGauss );
  }

  TH1D* etaSigmaHist = new TH1D( ("EtaSigma_"+timeTitle).c_str() , (timeTitle+";i#eta;Time Resolution #sigma (ns)").c_str() ,etaChannels, -85, 85);
  TH1D* etaMeanHist  = new TH1D( ("EtaMean_"+timeTitle ).c_str() , (timeTitle+";i#eta;Mean").c_str() ,etaChannels, -85, 85);
  TH1D* phiSigmaHist = new TH1D( ("PhiSigma_"+timeTitle).c_str() , (timeTitle+";i#phi;Time Resolution #sigma (ns)").c_str() ,phiChannels, 0, 360);
  TH1D* phiMeanHist  = new TH1D( ("PhiMean_"+timeTitle ).c_str() , (timeTitle+";i#phi;Mean").c_str() ,phiChannels, 0, 360);

  float etaSigMin = 999, etaMeanMin = 999, phiSigMin = 999, phiMeanMin = 999;
  for (int i=0; i<etaChannels; i++) {
    if (etaHists[i]->GetEntries() > 150 ) {
      etaSigmaHist->SetBinContent(i+1, etaSigma[i]);
      etaSigmaHist->SetBinError(i+1, etaSigEr[i]);
      etaMeanHist->SetBinContent(i+1, etaMean[i]);
      etaMeanHist->SetBinError(i+1, etaMeanEr[i]);
      if (etaSigma[i] - etaSigEr[i] < etaSigMin) etaSigMin = etaSigma[i] - etaSigEr[i];
      if (etaMean[i] - etaMeanEr[i] < etaMeanMin) etaMeanMin = etaMean[i] - etaMeanEr[i];
    }
  }
  for (int i=0; i<phiChannels; i++) {
    if (phiHists[i]->GetEntries() > 150 ) {
      phiSigmaHist->SetBinContent(i+1, phiSigma[i]);
      phiSigmaHist->SetBinError(i+1, phiSigEr[i]);
      phiMeanHist->SetBinContent(i+1, phiMean[i]);
      phiMeanHist->SetBinError(i+1, phiMeanEr[i]);
      if (phiSigma[i] - phiSigEr[i] < phiSigMin) phiSigMin = phiSigma[i] - phiSigEr[i];
      if (phiMean[i] - phiMeanEr[i] < phiMeanMin) phiMeanMin = phiMean[i] - phiMeanEr[i];
    }
  }

  etaSigmaHist->SetOption("ehist");
  etaMeanHist->SetOption("ehist");
  phiSigmaHist->SetOption("ehist");
  phiMeanHist->SetOption("ehist");

  etaSigmaHist->SetMinimum( etaSigMin );
  etaMeanHist->SetMinimum( etaMeanMin );
  phiSigmaHist->SetMinimum( phiSigMin );
  phiMeanHist->SetMinimum( phiMeanMin );

  SaveHist(etaSigmaHist);
  SaveHist(etaMeanHist);
  SaveHist(phiSigmaHist);
  SaveHist(phiMeanHist);
}



template <class U1>
U1 CutArray(U1 cuts[], int nSteps, U1 min, U1 max) {
  U1 stepSize = (max - min) / nSteps;
  for (int i=0; i<nSteps; i++) cuts[i] = min + i * stepSize;
  cuts[nSteps] = max+1; // Last element outside loop because StepSize might be rounded if int
  return stepSize; 
  // When looping through nSteps, should check quantity >= cuts[i] && quentity < cuts[i+1]. Note >= and <.
}


void hist2D( vector<vector<TH1F*> > hist, int zAxis, int xlen, float xmin, float xmax, int ylen, float ymin, float ymax, string histTitle, string Title, string xTitle, string yTitle) {
  // zAxis: 1 = Mean, 2 = Sigma, 0 = Entries
  TH2D *histFinal = new TH2D(histTitle.c_str(), string(Title+";"+xTitle+";"+yTitle).c_str() , xlen, xmin, xmax, ylen, ymin, ymax);

  if (zAxis == 0) {
    for (int i=0; i<xlen; i++) {
      for (int j=0; j<ylen; j++) {
        histFinal->SetBinContent(i+1, j+1, hist[i][j]->GetEntries());
      }
    }
  }


  else if (zAxis == 1) {
    TF1  *Fit[xlen][ylen];
    double param[xlen][ylen];
    double minParam = 5E7;
    bool singleGauss = 0;
    for (int i=0; i<xlen; i++) {
      for (int j=0; j<ylen; j++) {
        if ( hist[i][j]->GetEntries() > 150 ) {
          Fit[i][j] = Fitter( hist[i][j], singleGauss );
          param[i][j] = GetDGMean( Fit[i][j], singleGauss);
          if ( abs(param[i][j]) < 0.2 ) {
            histFinal->SetBinContent(i+1, j+1, param[i][j]);
            if ( param[i][j]<minParam ) minParam = param[i][j];
          }
          else histFinal->SetBinContent(i+1, j+1, -999);
        }
        else histFinal->SetBinContent(i+1, j+1, -999);
      }
    }
  histFinal->SetMinimum( minParam );
  }


  else if (zAxis == 2) {
    TF1  *Fit[xlen][ylen];
    double param[xlen][ylen];
    double minParam = 5E7;
    bool singleGauss = 0;
    for (int i=0; i<xlen; i++) {
      for (int j=0; j<ylen; j++) {
        if ( hist[i][j]->GetEntries() > 150 ) {
          Fit[i][j] = Fitter( hist[i][j], singleGauss );
          param[i][j] = GetDGSigma( Fit[i][j], singleGauss);
          if ( param[i][j] < 1.6 && param[i][j] > 0 ) {
            histFinal->SetBinContent(i+1, j+1, param[i][j]);
            if ( param[i][j]<minParam ) minParam = param[i][j];
          }
          else histFinal->SetBinContent(i+1, j+1, -999);
        }
        else histFinal->SetBinContent(i+1, j+1, -999);
      }
    }
  histFinal->SetMinimum( minParam );
  }


  histFinal->SetOption("colz");
  SaveHist(histFinal);
}



template <class T1>
void crystalTOF2(TH1F** hists, T1 cutArray[], int nBins, float energyCut_Abs, float energyUpperLimit, const vector<float> &energy, const vector<float> &eta, const vector<float> &phi, const vector<float> &transp, const vector<float> &time, int run, int xOpt, vector<int> &skip) {
  // Finds nearby, high-energy crystals, then fills histogram with their deltaT.
  // Returns a vector of the 2 indices, in decreasing energy order.
  if ( skip.empty() ) {
    auto maxEnergyItr = max_element( energy.begin(), energy.end() );
    float maxEnergy = *maxEnergyItr;
    if ( !(maxEnergy > energyCut_Abs) ) return;
    int crystal1Elt = distance(energy.begin(), maxEnergyItr);
    float energyCut_Rel = 0.7 * maxEnergy;
    float deltaRcut = 0.02;
    float deltaR1;
    vector<int> candidate_index; // Store indices of possible nearby crystals with 2nd highest energy
    for (unsigned int i=0; i<eta.size(); i++) {
      deltaR1 = sqrt( pow( eta[crystal1Elt] - eta[i] , 2)
                    + pow( phi[crystal1Elt] - phi[i] , 2) );
      if ( i!=crystal1Elt && deltaR1<deltaRcut && energy[i]>energyCut_Abs && energy[i]>energyCut_Rel ) candidate_index.push_back(i);
    }
    if ( candidate_index.empty() ) return;
    vector<int> candidate_energy;
    for (unsigned int i=0; i<candidate_index.size(); i++) candidate_energy.push_back( energy[candidate_index[i]] );
    int tempElt = distance( candidate_energy.begin(), max_element(candidate_energy.begin(), candidate_energy.end()) );
    int crystal2Elt = candidate_index[ tempElt ];
    skip.push_back(crystal1Elt);
    skip.push_back(crystal2Elt);
  }

  if (xOpt == 1) { // Energy
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(energy[skip[0]] > cutArray[i] && energy[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - time[skip[1]] );
      break;
    }
  }

  else if (xOpt == 2) { // Eta
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(eta[skip[0]] > cutArray[i] && eta[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - time[skip[1]] );
      break;
    }
  }

  else if (xOpt == 3) { // Phi
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(phi[skip[0]] > cutArray[i] && phi[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - time[skip[1]] );
      break;
    }
  }

  else if (xOpt == 4) { // Run
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(run > cutArray[i] && run < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - time[skip[1]] );
      break;
    }
  }

  else if (xOpt == 5) { // Transparency
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(transp[skip[0]] > cutArray[i] && transp[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - time[skip[1]] );
      break;
    }
  }
}




template <class T1>
void crystalTOF3(TH1F** hists, T1 cutArray[], int nBins, float energyCut_Abs, float energyUpperLimit, const vector<float> &energy, const vector<float> &eta, const vector<float> &phi, const vector<float> &transp, const vector<float> &time, int run, int xOpt, vector<int> &skip) {
  // Finds nearby, high-energy crystals, then fills histogram with their deltaT.
  // Returns a vector of the 3 indices, in decreasing energy order.
  if ( skip.empty() ) {
    auto maxEnergyItr = max_element( energy.begin(), energy.end() );
    float maxEnergy = *maxEnergyItr;
    if ( !(maxEnergy > energyCut_Abs) ) return;
    int crystal1Elt = distance(energy.begin(), maxEnergyItr);
    float energyCut_Rel = 0.7 * maxEnergy;
    float deltaRcut = 0.02;
    float deltaR1;
    vector<int> candidate_index; // Store indices of possible nearby crystals with 2nd highest energy
    for (unsigned int i=0; i<eta.size(); i++) {
      deltaR1 = sqrt( pow( eta[crystal1Elt] - eta[i] , 2)
                    + pow( phi[crystal1Elt] - phi[i] , 2) );
      if ( i!=crystal1Elt && deltaR1<deltaRcut && energy[i]>energyCut_Abs && energy[i]>energyCut_Rel ) candidate_index.push_back(i);
    }
    if ( candidate_index.empty() ) return;
    vector<int> candidate_energy;
    for (unsigned int i=0; i<candidate_index.size(); i++) candidate_energy.push_back( energy[candidate_index[i]] );
    int tempElt = distance( candidate_energy.begin(), max_element(candidate_energy.begin(), candidate_energy.end()) );
    int crystal2Elt = candidate_index[ tempElt ];
    // Repeat exact process for 3rd crystal
    candidate_index.clear();
    candidate_energy.clear();
    energyCut_Rel = 0.7 * energy[crystal2Elt];
    float deltaR2;
    for (unsigned int i=0; i<eta.size(); i++) {
      deltaR1 = sqrt( pow( eta[crystal1Elt] - eta[i] , 2)
                    + pow( phi[crystal1Elt] - phi[i] , 2) );
      deltaR2 = sqrt( pow( eta[crystal2Elt] - eta[i] , 2)
                    + pow( phi[crystal2Elt] - phi[i] , 2) );
      if ( i!=crystal1Elt && i!=crystal2Elt && (deltaR1<deltaRcut || deltaR2<deltaRcut) && energy[i]>energyCut_Abs && energy[i]>energyCut_Rel ) candidate_index.push_back(i);
    }
    if ( candidate_index.empty() ) return;
    for (unsigned int i=0; i<candidate_index.size(); i++) candidate_energy.push_back( energy[candidate_index[i]] );
    tempElt = distance( candidate_energy.begin(), max_element(candidate_energy.begin(), candidate_energy.end()) );
    int crystal3Elt = candidate_index[ tempElt ];
    skip.push_back(crystal1Elt);
    skip.push_back(crystal2Elt);
    skip.push_back(crystal3Elt);
  }
  // Combine times of crystals 2 and 3
  float newTime = (time[ skip[1] ]*energy[ skip[1] ] + time[ skip[2] ]*energy[ skip[2] ]) / (energy[ skip[2] ] + energy[ skip[3] ]);

  if (xOpt == 1) { // Energy
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(energy[skip[0]] > cutArray[i] && energy[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - newTime );
      break;
    }
  }

  else if (xOpt == 2) { // Eta
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(eta[skip[0]] > cutArray[i] && eta[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - newTime );
      break;
    }
  }

  else if (xOpt == 3) { // Phi
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(phi[skip[0]] > cutArray[i] && phi[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - newTime );
      break;
    }
  }

  else if (xOpt == 4) { // Run
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(run > cutArray[i] && run < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - newTime );
      break;
    }
  }

  else if (xOpt == 5) { // Transparency
    for (unsigned int i=0; i<nBins; i++) {
      if ( !(transp[skip[0]] > cutArray[i] && transp[skip[0]] < cutArray[i+1]) ) continue;
      hists[i]->Fill( time[skip[0]] - newTime );
      break;
    }
  }
}



TF1* AvgTimeGraph(TH1F** hists, int bins, float min, float max, string histTitle, string Title, string xTitle, string yTitle) {
  // hists is array of length bins
  TH1F *newHist = new TH1F(histTitle.c_str(), string(Title+";"+xTitle+";"+yTitle).c_str(), bins, min, max);
  for (int i=0; i<bins; i++) newHist->SetBinContent(i+1, hists[i]->GetMean() );
  return SaveHist(newHist, 1);
}




void Zee_DG() {
  string filename = "All2016.root";
  bool doCrystalTOF = 1, doRun = 1, doEventTime = 1, donPV = 1, doEtaPhi = 1, doTransparency = 1; // 0 to not re-generate plots

  TFile *inputfile = TFile::Open(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("ZeeTiming");

  float mass;
  unsigned int run;
  unsigned int nPV;
  unsigned int eventTime;
  vector<float> *ecalElectronRechit_EPtr = new vector<float>;
  vector<float> *ecalElectronRechit_EtaPtr = new vector<float>;
  vector<float> *ecalElectronRechit_PhiPtr = new vector<float>;
  vector<float> *ecalElectronRechit_transpCorrPtr = new vector<float>;
  vector<float> *ecalElectronRechit_calibT_legacyPtr = new vector<float>;
  vector<float> ecalElectronRechit_E; ecalElectronRechit_EPtr = &ecalElectronRechit_E;
  vector<float> ecalElectronRechit_Eta; ecalElectronRechit_EtaPtr = &ecalElectronRechit_Eta;
  vector<float> ecalElectronRechit_Phi; ecalElectronRechit_PhiPtr = &ecalElectronRechit_Phi;
  vector<float> ecalElectronRechit_transpCorr; ecalElectronRechit_transpCorrPtr = &ecalElectronRechit_transpCorr;
  vector<float> ecalElectronRechit_calibT_legacy; ecalElectronRechit_calibT_legacyPtr = &ecalElectronRechit_calibT_legacy;

  float t1;
  float t1_seed;
  float t1raw_seed;
  float t1calib_seed;
  float t1calib_seed_sept;
  int   ele1SeedIEta;
  int   ele1SeedIPhi;
  float ele1Eta;
  float ele1Phi;
  bool  ele1IsEB;
  float ele1Pt;
  float seed1_transpCorr;

  float t2;
  float t2_seed;
  float t2raw_seed;
  float t2calib_seed;
  float t2calib_seed_sept;
  int   ele2SeedIEta;
  int   ele2SeedIPhi;
  float ele2Eta;
  float ele2Phi;
  bool  ele2IsEB;
  float ele2Pt;
  float seed2_transpCorr;

  // Activate all branches to set addresses
  tree->SetBranchStatus("*",1);

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("nPV",&nPV);
  tree->SetBranchAddress("eventTime",&eventTime);
  tree->SetBranchAddress("mass",&mass);
  tree->SetBranchAddress("ecalElectronRechit_E",&ecalElectronRechit_EPtr);
  tree->SetBranchAddress("ecalElectronRechit_Eta",&ecalElectronRechit_EtaPtr);
  tree->SetBranchAddress("ecalElectronRechit_Phi",&ecalElectronRechit_PhiPtr);
  tree->SetBranchAddress("ecalElectronRechit_transpCorr",&ecalElectronRechit_transpCorrPtr);
  tree->SetBranchAddress("ecalElectronRechit_calibT_lagacy",&ecalElectronRechit_calibT_legacyPtr);

  tree->SetBranchAddress("t1",&t1);
  tree->SetBranchAddress("t1_seed",&t1_seed);
  tree->SetBranchAddress("t1raw_seed",&t1raw_seed);
  tree->SetBranchAddress("t1calib_seed",&t1calib_seed);
  tree->SetBranchAddress("t1calib_seed_sept",&t1calib_seed_sept);
  tree->SetBranchAddress("ele1SeedIEta",&ele1SeedIEta);
  tree->SetBranchAddress("ele1SeedIPhi",&ele1SeedIPhi);
  tree->SetBranchAddress("ele1Eta",&ele1Eta);
  tree->SetBranchAddress("ele1Phi",&ele1Phi);
  tree->SetBranchAddress("ele1IsEB",&ele1IsEB);
  tree->SetBranchAddress("ele1Pt",&ele1Pt);
  tree->SetBranchAddress("seed1_transpCorr",&seed1_transpCorr);

  tree->SetBranchAddress("t2",&t2);
  tree->SetBranchAddress("t2_seed",&t2_seed);
  tree->SetBranchAddress("t2raw_seed",&t2raw_seed);
  tree->SetBranchAddress("t2calib_seed",&t2calib_seed);
  tree->SetBranchAddress("t2calib_seed_sept",&t2calib_seed_sept);
  tree->SetBranchAddress("ele2SeedIEta",&ele2SeedIEta);
  tree->SetBranchAddress("ele2SeedIPhi",&ele2SeedIPhi);
  tree->SetBranchAddress("ele2Eta",&ele2Eta);
  tree->SetBranchAddress("ele2Phi",&ele2Phi);
  tree->SetBranchAddress("ele2IsEB",&ele2IsEB);
  tree->SetBranchAddress("ele2Pt",&ele2Pt);
  tree->SetBranchAddress("seed2_transpCorr",&seed2_transpCorr);

  // Activating the used branches only greatly speeds up looping through data
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("run",1);
  tree->SetBranchStatus("nPV",1);
  tree->SetBranchStatus("eventTime",1);

  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of Events in Sample: "<<nentries<<std::endl;

  // Find values of min and max
  float crystalEnergyMin = 10;
  float crystalEnergyMax = 105;
  float transpMin   = 0.4;
  float transpMax   = 3.4;
  unsigned int runMin = 10E7;
  unsigned int runMax = 0;
  unsigned int nPVMin = 1;
  unsigned int nPVMax = nPVMin;
  unsigned int eventTimeMin = 4E9;
  unsigned int eventTimeMax = 0;
  float tMin = -5;
  float tMax = 5;

  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    tree->GetEntry(iEntry);

    if      (run < runMin) runMin = run;
    else if (run > runMax) runMax = run; // the "else" should speed it up. If the data order is monotonic, remove it

    if      (eventTime < eventTimeMin) eventTimeMin = eventTime;
    else if (eventTime > eventTimeMax) eventTimeMax = eventTime;

    if (nPV > nPVMax) nPVMax = nPV;
  }

  // Construct (nSteps) intervals 
  int nSteps = 50;
  unsigned int runCuts[nSteps + 1];
  int stepSize = CutArray(runCuts, nSteps, runMin, runMax);

  unsigned int nPVCuts[nSteps + 1];
  int nPVStepSize = CutArray(nPVCuts, nSteps, nPVMin, nPVMax);

  unsigned int eventTimeCuts[nSteps + 1];
  int eventTimeStepSize = CutArray(eventTimeCuts, nSteps, eventTimeMin, eventTimeMax);

  float transpCuts[nSteps+1];
  float transpStepSize = CutArray(transpCuts, nSteps, transpMin, transpMax);

  float tCuts[nSteps + 1];
  float tStepSize = CutArray(tCuts, nSteps, tMin, tMax);

  float crystalEnergyCuts[nSteps + 1];
  float crystalEnergyStepSize = CutArray(crystalEnergyCuts, nSteps, crystalEnergyMin, crystalEnergyMax);

  // Declare Histograms
  TH1F *histRun_t[nSteps];
  TH1F *histRun_tseed[nSteps];
  TH1F *histRun_trawseed[nSteps];
  TH1F *histRun_tcalibseed[nSteps];
  TH1F *histRun_tcalibseedsept[nSteps];

  TH1F *histnPV_t[nSteps];

  TH1F *histEventTime_t[nSteps];

  TH1F *histTransp1_t[nSteps];
  TH1F *histTransp2_t[nSteps];
  TH1F *histTransp1_tseed[nSteps];
  TH1F *histTransp2_tseed[nSteps];
  TH1F *histTransp1_trawseed[nSteps];
  TH1F *histTransp2_trawseed[nSteps];
  TH1F *histTransp1_tcalibseed[nSteps];
  TH1F *histTransp2_tcalibseed[nSteps];
  TH1F *histTransp1_tcalibseedsept[nSteps];
  TH1F *histTransp2_tcalibseedsept[nSteps];

  // Below are 2D vector arrays of TH1F histograms
  vector< vector<TH1F*> > histTransp1Transp2_t (nSteps, vector<TH1F*>(nSteps)); // transp2 vs transp 1 (sig & mean on Z)
  vector< vector<TH1F*> > histTransp1Transp2_tseed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1Transp2_trawseed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1Transp2_tcalibseed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1Transp2_tcalibseedsept (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1t1 (nSteps, vector<TH1F*>(nSteps)); // t1 vs Transp1 (nEvents on Z)
  vector< vector<TH1F*> > histTransp2t2 (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1t1seed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp2t2seed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1t1rawseed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp2t2rawseed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1t1calibseed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp2t2calibseed (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp1t1calibseedsept (nSteps, vector<TH1F*>(nSteps));
  vector< vector<TH1F*> > histTransp2t2calibseedsept (nSteps, vector<TH1F*>(nSteps));
  TH1F *histTransp1_t1[nSteps]; // average t1 value for each bin of Transp1
  TH1F *histTransp2_t2[nSteps];
  TH1F *histTransp1_t1seed[nSteps];
  TH1F *histTransp2_t2seed[nSteps];
  TH1F *histTransp1_t1rawseed[nSteps];
  TH1F *histTransp2_t2rawseed[nSteps];
  TH1F *histTransp1_t1calibseed[nSteps];
  TH1F *histTransp2_t2calibseed[nSteps];
  TH1F *histTransp1_t1calibseedsept[nSteps];
  TH1F *histTransp2_t2calibseedsept[nSteps];

  int BinSize = 10; // 1, 2, 5, or 10
  int phiChannels = 360;
  int etaChannels = 170;
  int phiBins = phiChannels/BinSize;
  int etaBins = etaChannels/BinSize;
  float etaCuts[etaChannels + 1];
  float etaStepSize = CutArray(etaCuts, etaChannels, float(-1.479), float(1.479));
  float phiCuts[phiChannels + 1];
  float phiStepSize = CutArray(phiCuts, phiChannels, float(-3.14), float(3.14));
  vector< vector<TH1F*> > histEtaPhi_t (phiBins, vector<TH1F*>(etaBins)); // Phi from 0 to 360 inclusive. Eta from -85 to 85 inclusive. Bins of 3.
  TH1F *histEta_t[etaChannels];
  TH1F *histPhi_t[phiChannels];
  TH1F *histEta_tseed[etaChannels];
  TH1F *histPhi_tseed[phiChannels];
  TH1F *histEta_trawseed[etaChannels];
  TH1F *histPhi_trawseed[phiChannels];
  TH1F *histEta_tcalibseed[etaChannels];
  TH1F *histPhi_tcalibseed[phiChannels];
  TH1F *histEta_tcalibseedsept[etaChannels];
  TH1F *histPhi_tcalibseedsept[phiChannels];

  TH1F *histCrystalTOFEnergy[nSteps];
  TH1F *histCrystalTOFTransp[nSteps];
  TH1F *histCrystalTOFRun[nSteps];
  TH1F *histCrystalTOFPhi[phiChannels];
  TH1F *histCrystalTOFEta[etaChannels];

  // Initialize Histograms
  for(int i=0; i<nSteps; i++) {
    histCrystalTOFEnergy[i] = new TH1F( Form("histCrystalTOFEnergy[%d]",i),";t_{Crystal 1}-t_{Crystal 2};Entries", 120, -1, 1);
    histCrystalTOFTransp[i] = new TH1F( Form("histCrystalTOFTransp[%d]",i),";t_{Crystal 1}-t_{Crystal 2};Entries", 120, -1, 1);
    histCrystalTOFRun[i] = new TH1F( Form("histCrystalTOFRun[%d]",i),";t_{Crystal 1}-t_{Crystal 2};Entries", 120, -1, 1);

    histRun_t[i]     = new TH1F( Form("histRun_t[%d]",i) ,";t_{1}-t_{2};Entries", 120, -3, 3);
    histRun_tseed[i] = new TH1F( Form("histRun_tseed[%d]",i) ,";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histRun_trawseed[i] = new TH1F( Form("histRun_trawseed[%d]",i) ,";t_{1}-t_{2} raw seed;Entries", 120, -5, 5);
    histRun_tcalibseed[i] = new TH1F( Form("histRun_tcalibseed[%d]",i) ,";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histRun_tcalibseedsept[i] = new TH1F( Form("histRun_tcalibseedsept[%d]",i) ,";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);

    histnPV_t[i] = new TH1F( Form("histnPV_t[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);

    histEventTime_t[i] = new TH1F( Form("histEventTime[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);

    histTransp1_t[i] = new TH1F( Form("histTransp1_t[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);
    histTransp2_t[i] = new TH1F( Form("histTransp2_t[%d]",i), ";t_{1}-t_{2};Entries", 120, -3, 3);
    histTransp1_tseed[i] = new TH1F( Form("histTransp1_tseed[%d]",i), ";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histTransp2_tseed[i] = new TH1F( Form("histTransp2_tseed[%d]",i), ";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histTransp1_trawseed[i] = new TH1F( Form("histTransp1_trawseed[%d]",i), ";t_{1}-t_{2} raw seed;Entries", 120, -5, 5);
    histTransp2_trawseed[i] = new TH1F( Form("histTransp2_trawseed[%d]",i), ";t_{1}-t_{2} raw seed;Entries", 120, -5, 5);
    histTransp1_tcalibseed[i] = new TH1F( Form("histTransp1_tcalibseed[%d]",i), ";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histTransp2_tcalibseed[i] = new TH1F( Form("histTransp2_tcalibseed[%d]",i), ";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histTransp1_tcalibseedsept[i] = new TH1F( Form("histTransp1_tcalibseedsept[%d]",i), ";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);
    histTransp2_tcalibseedsept[i] = new TH1F( Form("histTransp2_tcalibseedsept[%d]",i), ";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);

    for (int j=0; j<nSteps; j++) {
      histTransp1Transp2_t[i][j] = new TH1F( Form("histTransp1Transp2_t[%d][%d]",i,j), ";t_{1}-t_{2};Entries", 120, -3, 3);
      histTransp1Transp2_tseed[i][j] = new TH1F( Form("histTransp1Transp2_tseed[%d][%d]",i,j), ";t_{1}-t_{2} Seed;Entries", 120, -3, 3);
      histTransp1Transp2_trawseed[i][j] = new TH1F( Form("histTransp1Transp2_trawseed[%d][%d]",i,j), ";t_{1}-t_{2} Raw Seed;Entries", 120, -5, 5);
      histTransp1Transp2_tcalibseed[i][j] = new TH1F( Form("histTransp1Transp2_tcalibseed[%d][%d]",i,j), ";t_{1}-t_{2} Calib Seed;Entries", 120, -3, 3);
      histTransp1Transp2_tcalibseedsept[i][j] = new TH1F( Form("histTransp1Transp2_tcalibseedsept[%d][%d]",i,j), ";t_{1}-t_{2} Calib Seed Sept;Entries", 120, -3, 3);

      histTransp1t1[i][j] = new TH1F( Form("histTransp1t1[%d][%d]",i,j), ";t_{1};Entries", 120, -3, 3);
      histTransp2t2[i][j] = new TH1F( Form("histTransp2t2[%d][%d]",i,j), ";t_{2};Entries", 120, -3, 3);;
      histTransp1t1seed[i][j] = new TH1F( Form("histTransp1t1seed[%d][%d]",i,j), ";t_{1} Seed;Entries", 120, -3, 3);
      histTransp2t2seed[i][j] = new TH1F( Form("histTransp2t2seed[%d][%d]",i,j), ";t_{2} Seed;Entries", 120, -3, 3);;
      histTransp1t1rawseed[i][j] = new TH1F( Form("histTransp1t1rawseed[%d][%d]",i,j), ";t_{1} Raw Seed;Entries", 120, -5, 5);
      histTransp2t2rawseed[i][j] = new TH1F( Form("histTransp2t2rawseed[%d][%d]",i,j), ";t_{2} Raw Seed;Entries", 120, -5, 5);;
      histTransp1t1calibseed[i][j] = new TH1F( Form("histTransp1t1calibseed[%d][%d]",i,j), ";t_{1} Calib Seed;Entries", 120, -3, 3);
      histTransp2t2calibseed[i][j] = new TH1F( Form("histTransp2t2calibseed[%d][%d]",i,j), ";t_{2} Calib Seed;Entries", 120, -3, 3);;
      histTransp1t1calibseedsept[i][j] = new TH1F( Form("histTransp1t1calibseedsept[%d][%d]",i,j), ";t_{1} Calib Seed Sept;Entries", 120, -3, 3);
      histTransp2t2calibseedsept[i][j] = new TH1F( Form("histTransp2t2calibseedsept[%d][%d]",i,j), ";t_{2} Calib Seed Sept;Entries", 120, -3, 3);;
    }

      histTransp1_t1[i] = new TH1F( Form("histTransp1_t1[%d]",i), ";t_{1};Entries", 120, -3, 3); 
      histTransp2_t2[i] = new TH1F( Form("histTransp2_t2[%d]",i), ";t_{2};Entries", 120, -3, 3);
      histTransp1_t1seed[i] = new TH1F( Form("histTransp1_t1seed[%d]",i), ";t_{1} Seed;Entries", 12000, -3, 3);
      histTransp2_t2seed[i] = new TH1F( Form("histTransp2_t2seed[%d]",i), ";t_{2} Seed;Entries", 120, -3, 3);
      histTransp1_t1rawseed[i] = new TH1F( Form("histTransp1_t1rawseed[%d]",i), ";t_{1} Raw Seed;Entries", 120, -5, 5);
      histTransp2_t2rawseed[i] = new TH1F( Form("histTransp2_t2rawseed[%d]",i), ";t_{2} Raw Seed;Entries", 120, -5, 5);
      histTransp1_t1calibseed[i] = new TH1F( Form("histTransp1_t1calibseed[%d]",i), ";t_{1} Calib Seed;Entries", 120, -3, 3);
      histTransp2_t2calibseed[i] = new TH1F( Form("histTransp2_t2calibseed[%d]",i), ";t_{2} Calib Seed;Entries", 120, -3, 3);
      histTransp1_t1calibseedsept[i] = new TH1F( Form("histTransp1_t1calibseedsept[%d]",i), ";t_{1} Calib Seed Sept;Entries", 120, -3, 3);
      histTransp2_t2calibseedsept[i] = new TH1F( Form("histTransp2_t2calibseedsept[%d]",i), ";t_{2} Calib Seed Sept;Entries", 120, -3, 3);
  }
  for (int i=0; i<phiBins; i++){
    for (int j=0; j<etaBins; j++){
      histEtaPhi_t[i][j] = new TH1F( Form("histEtaPhi_t[%d][%d]",i,j) ,";t_{1}-t_{2};Entries", 120, -3, 3);
    }
  }
  for (int i=0; i<etaChannels; i++) {
    histEta_t[i] = new TH1F( Form("histEta_t[%d]",i),";t_{1}-t_{2};Entries", 120, -3, 3);
    histEta_tseed[i] = new TH1F( Form("histEta_tseed[%d]",i),";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histEta_trawseed[i] = new TH1F( Form("histEta_trawseed[%d]",i),";t_{1}-t_{2} raw seed;Entries", 120, -5, 5);
    histEta_tcalibseed[i] = new TH1F( Form("histEtacalibseed_t[%d]",i),";t_{1}-t_{2} calib seed;Entries", 120, -3, 3);
    histEta_tcalibseedsept[i] = new TH1F( Form("histEta_tcalibseedsept[%d]",i),";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);

    histCrystalTOFEta[i] = new TH1F( Form("histCrystalTOFEta[%d]",i),";t_{Crystal 1}-t_{Crystal 2};Entries", 120, -1, 1);
  }
  for (int i=0; i<phiChannels; i++) {
    histPhi_t[i] = new TH1F( Form("histPhi_t[%d]",i),";t_{1}-t_{2};Entries", 120, -3, 3);
    histPhi_tseed[i] = new TH1F( Form("histPhi_tseed[%d]",i),";t_{1}-t_{2} seed;Entries", 120, -3, 3);
    histPhi_trawseed[i] = new TH1F( Form("histPhi_trawseed[%d]",i),";t_{1}-t_{2} raw seed;Entries", 120, -5, 5);
    histPhi_tcalibseed[i] = new TH1F( Form("histPhi_tcalibseed[%d]",i),";t_{1}-t_{2 calib seed};Entries", 120, -3, 3);
    histPhi_tcalibseedsept[i] = new TH1F( Form("histPhi_tcalibseedsept[%d]",i),";t_{1}-t_{2} calib seed sept;Entries", 120, -3, 3);

    histCrystalTOFPhi[i] = new TH1F( Form("histCrystalTOFPhi[%d]",i),";t_{Crystal 1}-t_{Crystal 2};Entries", 120, -1, 1);
  }


  TFile *file = TFile::Open(("output"+filename).c_str(), "RECREATE");
  file->cd();

  // FOR EACH PROCESS:
  // if statement to check if process is enabled
  // activate necessary branches for process
  // commence event loop and fill histograms
  // exit loop and create plots

  

  cout<<"\nCrystal TOF"<<endl;
  if (doCrystalTOF) {
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("mass",1);
    tree->SetBranchStatus("ecalElectronRechit_E",1);
    tree->SetBranchStatus("ecalElectronRechit_Eta",1);
    tree->SetBranchStatus("ecalElectronRechit_Phi",1);
    tree->SetBranchStatus("ecalElectronRechit_transpCorr",1);
    tree->SetBranchStatus("ecalElectronRechit_calibT_lagacy",1);

    tree->SetBranchStatus("ele1IsEB",1);
    tree->SetBranchStatus("ele1Pt",1);

    tree->SetBranchStatus("ele2IsEB",1);
    tree->SetBranchStatus("ele2Pt",1);

    for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
      if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);

      if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything

      vector<int> energyElts;
      // change crystalTOF3 to crystalTOF2 to do the 2-crystal TOF analysis
      crystalTOF2(histCrystalTOFEnergy, crystalEnergyCuts, nSteps, crystalEnergyMin, crystalEnergyMax, ecalElectronRechit_E, ecalElectronRechit_Eta, ecalElectronRechit_Phi, ecalElectronRechit_transpCorr, ecalElectronRechit_calibT_legacy, run, 1, energyElts );
      crystalTOF2(histCrystalTOFEta, etaCuts, etaChannels, crystalEnergyMin, crystalEnergyMax, ecalElectronRechit_E, ecalElectronRechit_Eta, ecalElectronRechit_Phi, ecalElectronRechit_transpCorr, ecalElectronRechit_calibT_legacy, run, 2, energyElts);
      crystalTOF2(histCrystalTOFPhi, phiCuts, phiChannels, crystalEnergyMin, crystalEnergyMax, ecalElectronRechit_E, ecalElectronRechit_Eta, ecalElectronRechit_Phi, ecalElectronRechit_transpCorr, ecalElectronRechit_calibT_legacy, run, 3, energyElts);
      crystalTOF2(histCrystalTOFRun, runCuts, nSteps, crystalEnergyMin, crystalEnergyMax, ecalElectronRechit_E, ecalElectronRechit_Eta, ecalElectronRechit_Phi, ecalElectronRechit_transpCorr, ecalElectronRechit_calibT_legacy, run, 4, energyElts);
      crystalTOF2(histCrystalTOFTransp, transpCuts, nSteps, crystalEnergyMin, crystalEnergyMax, ecalElectronRechit_E, ecalElectronRechit_Eta, ecalElectronRechit_Phi, ecalElectronRechit_transpCorr, ecalElectronRechit_calibT_legacy, run, 5, energyElts);
    }

    cout << "Generating Crystal TOF Plots..." << endl;
    SigMeanTGraph(file, histCrystalTOFEnergy, "histCrystalTOFEnergy", crystalEnergyCuts, crystalEnergyStepSize, nSteps, "", "Highest Crystal Energy (GeV)", "crystalTOFEnergy");
    SigMeanTGraph(file, histCrystalTOFEta, "histCrystalTOFEta", etaCuts, etaStepSize, etaChannels, "", "#eta", "crystalTOFEta");
    SigMeanTGraph(file, histCrystalTOFPhi, "histCrystalTOFPhi", phiCuts, phiStepSize, phiChannels, "", "#phi", "crystalTOFPhi");
    SigMeanTGraph(file, histCrystalTOFRun, "histCrystalTOFRun", runCuts, stepSize, nSteps, "", "Run Number", "crystalTOFRun");
    SigMeanTGraph(file, histCrystalTOFTransp, "histCrystalTOFTransp", transpCuts, transpStepSize, nSteps, "", "Transparency", "crystalTOFTransp");
    TH1F *histCrystalTOF = new TH1F( "histCrystalTOF",";t_{Crystal 1}-t_{Crystal 2};Entries", 120, -1, 1);
    for (unsigned int i=0; i<nSteps; i++) histCrystalTOF->Add(histCrystalTOFEnergy[i]);
    file->WriteTObject(histCrystalTOF, "histCrystalTOF","WriteDelete");
  }



  cout<<"\nRun"<<endl;
  if (doRun) {
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("mass",1);

    tree->SetBranchStatus("t1",1);
    tree->SetBranchStatus("t1_seed",1);
    tree->SetBranchStatus("t1raw_seed",1);
    tree->SetBranchStatus("t1calib_seed",1);
    tree->SetBranchStatus("t1calib_seed_sept",1);
    tree->SetBranchStatus("ele1IsEB",1);
    tree->SetBranchStatus("ele1Pt",1);

    tree->SetBranchStatus("t2",1);
    tree->SetBranchStatus("t2_seed",1);
    tree->SetBranchStatus("t2raw_seed",1);
    tree->SetBranchStatus("t2calib_seed",1);
    tree->SetBranchStatus("t2calib_seed_sept",1);
    tree->SetBranchStatus("ele2IsEB",1);
    tree->SetBranchStatus("ele2Pt",1);

    for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
      if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);

      if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything

      for (int i=0; i<nSteps; i++) {
        if( !(run>=runCuts[i] && run<runCuts[i+1]) ) continue;
        histRun_t[i]->Fill( t1-t2 ); 
        histRun_tseed[i]->Fill( t1_seed-t2_seed );
        histRun_trawseed[i]->Fill( t1raw_seed-t2raw_seed );
        histRun_tcalibseed[i]->Fill( t1calib_seed-t2calib_seed );
        histRun_tcalibseedsept[i]->Fill( t1calib_seed_sept-t2calib_seed_sept );
        break; // Don't bother to check other cases in the for loop
      }
    }

    std::cout<<"Generating Run Plots..."<<std::endl;
    SigMeanTGraph(file, histRun_t, "histRun_t", runCuts, stepSize, nSteps, "", "Run Number", "run_t");
    SigMeanTGraph(file, histRun_tseed, "histRun_tseed", runCuts, stepSize, nSteps, "Seed", "Run Number", "run_tseed");
    SigMeanTGraph(file, histRun_trawseed, "histRun_trawseed", runCuts, stepSize, nSteps, "Raw Seed", "Run Number", "run_trawseed");
    SigMeanTGraph(file, histRun_tcalibseed, "histRun_tcalibseed", runCuts, stepSize, nSteps, "Calib Seed", "Run Number", "run_tcalibseed");
    SigMeanTGraph(file, histRun_tcalibseedsept, "histRun_tcalibseedsept", runCuts, stepSize, nSteps, "Calib Seed Sept", "Run Number", "run_tcalibseedsept");
  }



  cout<<"\nnPV"<<endl;
  if (donPV) {
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("nPV",1);
    tree->SetBranchStatus("mass",1);

    tree->SetBranchStatus("t1",1);
    tree->SetBranchStatus("ele1IsEB",1);
    tree->SetBranchStatus("ele1Pt",1);

    tree->SetBranchStatus("t2",1);
    tree->SetBranchStatus("ele2IsEB",1);
    tree->SetBranchStatus("ele2Pt",1);

    for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
      if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);

      if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything

      for (int i=0; i<nSteps; i++) {
        if( !(nPV>=nPVCuts[i] && nPV<nPVCuts[i+1]) ) continue;
          histnPV_t[i]->Fill( t1-t2 ); 
          break;
      }
    }
    std::cout<<"Generating nPV Plots..."<<std::endl;
    SigMeanTGraph(file, histnPV_t, "histnPV_t", nPVCuts, nPVStepSize, nSteps, "", "nPV", "nPV_t");
  }



  cout<<"\neventTime"<<endl;
  if (doEventTime) {
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("eventTime",1);
    tree->SetBranchStatus("mass",1);

    tree->SetBranchStatus("t1",1);
    tree->SetBranchStatus("ele1IsEB",1);
    tree->SetBranchStatus("ele1Pt",1);

    tree->SetBranchStatus("t2",1);
    tree->SetBranchStatus("ele2IsEB",1);
    tree->SetBranchStatus("ele2Pt",1);

    for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
      if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);

      if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything

      for (int i=0; i<nSteps; i++) {
        if( !(eventTime>=eventTimeCuts[i] && eventTime<eventTimeCuts[i+1]) ) continue;
          histEventTime_t[i]->Fill( t1-t2 );
          break;
      }
    }
    std::cout<<"Generating EventTime Plots..."<<std::endl;
    SigMeanTGraph(file, histEventTime_t, "histEventTime_t", eventTimeCuts, eventTimeStepSize, nSteps, "", "Event Time", "eventTime_t");
  }



  cout<<"\nEta & Phi"<<endl;
  if (doEtaPhi) {
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("mass",1);

    tree->SetBranchStatus("t1",1);
    tree->SetBranchStatus("t1_seed",1);
    tree->SetBranchStatus("t1raw_seed",1);
    tree->SetBranchStatus("t1calib_seed",1);
    tree->SetBranchStatus("t1calib_seed_sept",1);
    tree->SetBranchStatus("ele1SeedIEta",1);
    tree->SetBranchStatus("ele1SeedIPhi",1);
    tree->SetBranchStatus("ele1Eta",1);
    tree->SetBranchStatus("ele1Phi",1);
    tree->SetBranchStatus("ele1IsEB",1);
    tree->SetBranchStatus("ele1Pt",1);

    tree->SetBranchStatus("t2",1);
    tree->SetBranchStatus("t2_seed",1);
    tree->SetBranchStatus("t2raw_seed",1);
    tree->SetBranchStatus("t2calib_seed",1);
    tree->SetBranchStatus("t2calib_seed_sept",1);
    tree->SetBranchStatus("ele2SeedIEta",1);
    tree->SetBranchStatus("ele2SeedIPhi",1);
    tree->SetBranchStatus("ele2Eta",1);
    tree->SetBranchStatus("ele2Phi",1);
    tree->SetBranchStatus("ele2IsEB",1);
    tree->SetBranchStatus("ele2Pt",1);

    for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
      if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);

      if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything

      if( abs(ele1Eta)<1.47 && abs(ele2Eta)<1.47 ) { // Electrons in barrel
        if( ele1SeedIEta > 0) ele1SeedIEta -= 1; // there are no events for ele1SeedIEta == 0 so we shift >0 entries down
        histEtaPhi_t[ (ele1SeedIPhi-1)/BinSize ][ (ele1SeedIEta+85)/BinSize ]->Fill(t1-t2);
        histEta_t[ (ele1SeedIEta+85) ]->Fill( t1-t2 );
        histPhi_t[ (ele1SeedIPhi-1)  ]->Fill( t1-t2 );
        histEta_tseed[ (ele1SeedIEta+85) ]->Fill( t1_seed-t2_seed );
        histPhi_tseed[ (ele1SeedIPhi-1)  ]->Fill( t1_seed-t2_seed );
        histEta_trawseed[ (ele1SeedIEta+85) ]->Fill( t1raw_seed-t2raw_seed );
        histPhi_trawseed[ (ele1SeedIPhi-1)  ]->Fill( t1raw_seed-t2raw_seed );
        histEta_tcalibseed[ (ele1SeedIEta+85) ]->Fill( t1calib_seed-t2calib_seed );
        histPhi_tcalibseed[ (ele1SeedIPhi-1)  ]->Fill( t1calib_seed-t2calib_seed );
        histEta_tcalibseedsept[ (ele1SeedIEta+85) ]->Fill( t1calib_seed_sept-t2calib_seed_sept );
        histPhi_tcalibseedsept[ (ele1SeedIPhi-1)  ]->Fill( t1calib_seed_sept-t2calib_seed_sept );
        if( ele1SeedIEta >= 0) ele1SeedIEta += 1; // Reset value when leaving loop
      }
    }

    std::cout<<"Generating 1D Eta and Phi Plots..."<<std::endl;
    EtaPhi1D( etaChannels, phiChannels, histEta_t,              histPhi_t,              "t" );
    EtaPhi1D( etaChannels, phiChannels, histEta_tseed,          histPhi_tseed,          "tseed" );
    EtaPhi1D( etaChannels, phiChannels, histEta_trawseed,       histPhi_trawseed,       "trawseed" );
    EtaPhi1D( etaChannels, phiChannels, histEta_tcalibseed,     histPhi_tcalibseed,     "tcalibseed" );
    EtaPhi1D( etaChannels, phiChannels, histEta_tcalibseedsept, histPhi_tcalibseedsept, "tcalibseedsept" );

    cout<<"Generating mean:eta:phi Plot..."<<endl;
    hist2D( histEtaPhi_t, 1, phiBins,0,360, etaBins,-85,85, "EtaPhiMean_t", "Mean", "i#phi", "i#eta" );

    cout<<"Generating sigma:eta:phi Plot..."<<endl;
    hist2D( histEtaPhi_t, 2, phiBins,0,360, etaBins,-85,85, "EtaPhiSigma_t","#sigma", "i#phi", "i#eta" );
  }



  cout<<"\nTransparency"<<endl;
  if (doTransparency) {
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("mass",1);

    tree->SetBranchStatus("t1",1);
    tree->SetBranchStatus("t1_seed",1);
    tree->SetBranchStatus("t1raw_seed",1);
    tree->SetBranchStatus("t1calib_seed",1);
    tree->SetBranchStatus("t1calib_seed_sept",1);
    tree->SetBranchStatus("ele1IsEB",1);
    tree->SetBranchStatus("ele1Pt",1);
    tree->SetBranchStatus("seed1_transpCorr",1);

    tree->SetBranchStatus("t2",1);
    tree->SetBranchStatus("t2_seed",1);
    tree->SetBranchStatus("t2raw_seed",1);
    tree->SetBranchStatus("t2calib_seed",1);
    tree->SetBranchStatus("t2calib_seed_sept",1);
    tree->SetBranchStatus("ele2IsEB",1);
    tree->SetBranchStatus("ele2Pt",1);
    tree->SetBranchStatus("seed2_transpCorr",1);

    for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
      if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
      tree->GetEntry(iEntry);

      if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything

      for (int i=0; i<nSteps; i++) {
        if( !(seed1_transpCorr>=transpCuts[i] && seed1_transpCorr<transpCuts[i+1]) ) continue;
        histTransp1_t[i]->Fill( t1-t2 );
        histTransp1_tseed[i]->Fill( t1_seed-t2_seed );
        histTransp1_trawseed[i]->Fill( t1raw_seed-t2raw_seed );
        histTransp1_tcalibseed[i]->Fill( t1calib_seed-t2calib_seed );
        histTransp1_tcalibseedsept[i]->Fill( t1calib_seed_sept-t2calib_seed_sept );
        
        histTransp1_t1[i]->Fill( t1 );
        histTransp1_t1seed[i]->Fill( t1_seed );
        histTransp1_t1rawseed[i]->Fill( t1raw_seed );
        histTransp1_t1calibseed[i]->Fill( t1calib_seed );
        histTransp1_t1calibseedsept[i]->Fill( t1calib_seed_sept );

        break; // Don't bother to check other cases in the for loop
      }
      for (int i=0; i<nSteps; i++) {
        if( !(seed2_transpCorr>=transpCuts[i] && seed2_transpCorr<transpCuts[i+1]) ) continue;
        histTransp2_t[i]->Fill( t1-t2 );
        histTransp2_tseed[i]->Fill( t1_seed-t2_seed );
        histTransp2_trawseed[i]->Fill( t1raw_seed-t2raw_seed );
        histTransp2_tcalibseed[i]->Fill( t1calib_seed-t2calib_seed );
        histTransp2_tcalibseedsept[i]->Fill( t1calib_seed_sept-t2calib_seed_sept );

        histTransp2_t2[i]->Fill( t2 );
        histTransp2_t2seed[i]->Fill( t2_seed );
        histTransp2_t2rawseed[i]->Fill( t2raw_seed );
        histTransp2_t2calibseed[i]->Fill( t2calib_seed );
        histTransp2_t2calibseedsept[i]->Fill( t2calib_seed_sept );

        break; // Don't bother to check other cases in the for loop
      }
      int brk = 0;
      for (int i=0; i<nSteps; i++) {
        for (int j=0; j<nSteps; j++) {
          if( seed1_transpCorr>=transpCuts[i] && seed1_transpCorr<transpCuts[i+1] &&
              seed2_transpCorr>=transpCuts[j] && seed2_transpCorr<transpCuts[j+1]) {
            histTransp1Transp2_t[i][j]->Fill( t1-t2 );
            histTransp1Transp2_tseed[i][j]->Fill( t1_seed-t2_seed );
            histTransp1Transp2_trawseed[i][j]->Fill( t1raw_seed-t2raw_seed );
            histTransp1Transp2_tcalibseed[i][j]->Fill( t1calib_seed-t2calib_seed );
            histTransp1Transp2_tcalibseedsept[i][j]->Fill( t1calib_seed_sept-t2calib_seed_sept );
            brk += 1;
          }
          if( seed1_transpCorr>=transpCuts[i] && seed1_transpCorr<transpCuts[i+1] && t1>=tCuts[j] && t1<tCuts[j+1]) {
            histTransp1t1[i][j]->Fill( t1 );// Fill content doesn't matter, just want nEntries
            histTransp1t1seed[i][j]->Fill( t1_seed );
            histTransp1t1rawseed[i][j]->Fill( t1raw_seed );
            histTransp1t1calibseed[i][j]->Fill( t1calib_seed );
            histTransp1t1calibseedsept[i][j]->Fill( t1calib_seed_sept );
            brk += 1;
          }
          if( seed2_transpCorr>=transpCuts[i] && seed2_transpCorr<transpCuts[i+1] && t2>=tCuts[j] && t2<tCuts[j+1]) {
            histTransp2t2[i][j]->Fill( t2 );// Fill content doesn't matter, just want nEntries
            histTransp2t2seed[i][j]->Fill( t2_seed );
            histTransp2t2rawseed[i][j]->Fill( t2raw_seed );
            histTransp2t2calibseed[i][j]->Fill( t2calib_seed );
            histTransp2t2calibseedsept[i][j]->Fill( t2calib_seed_sept );
            brk += 1;
          }
        }
        if (brk == 3) break;
      }
    }
    std::cout<<"Generating Transparency Plots..."<<std::endl;
    SigMeanTGraph(file, histTransp1_t, "histTransp1_t", transpCuts, transpStepSize, nSteps, "", "Seed 1 Transparency", "transp1_t");
    SigMeanTGraph(file, histTransp2_t, "histTransp2_t", transpCuts, transpStepSize, nSteps, "", "Seed 2 Transparency", "transp2_t");
    SigMeanTGraph(file, histTransp1_tseed, "histTransp1_tseed", transpCuts, transpStepSize, nSteps, "Seed", "Seed 1 Transparency", "transp1_tseed");
    SigMeanTGraph(file, histTransp2_tseed, "histTransp2_tseed", transpCuts, transpStepSize, nSteps, "Seed", "Seed 2 Transparency", "transp2_tseed");
    SigMeanTGraph(file, histTransp1_trawseed, "histTransp1_trawseed", transpCuts, transpStepSize, nSteps, "Raw Seed", "Seed 1 Transparency", "transp1_trawseed");
    SigMeanTGraph(file, histTransp2_trawseed, "histTransp2_trawseed", transpCuts, transpStepSize, nSteps, "Raw Seed", "Seed 2 Transparency", "transp2_trawseed");
    SigMeanTGraph(file, histTransp1_tcalibseed, "histTransp1_tcalibseed", transpCuts, transpStepSize, nSteps, "Calib Seed", "Seed 1 Transparency", "transp1_tcalibseed");
    SigMeanTGraph(file, histTransp2_tcalibseed, "histTransp2_tcalibseed", transpCuts, transpStepSize, nSteps, "Calib Seed", "Seed 2 Transparency", "transp2_tcalibseed");
    SigMeanTGraph(file, histTransp1_tcalibseedsept, "histTransp1_tcalibseedsept", transpCuts, transpStepSize, nSteps, "Calib Seed Sept", "Seed 1 Transparency", "transp1_tcalibseedsept");
    SigMeanTGraph(file, histTransp2_tcalibseedsept, "histTransp2_tcalibseedsept", transpCuts, transpStepSize, nSteps, "Calib Seed Sept", "Seed 2 Transparency", "transp2_tcalibseedsept");


    cout<<"Generating Avg Time vs Transparency Plots..."<<endl;
    TF1 *fit1_t = AvgTimeGraph(histTransp1_t1, nSteps, transpMin, transpMax, "Transp1_t1", "", "Seed 1 Transparency", "Average t_{1}");
    TF1 *fit2_t = AvgTimeGraph(histTransp2_t2, nSteps, transpMin, transpMax, "Transp2_t2", "", "Seed 2 Transparency", "Average t_{2}");
    TF1 *fit1_tseed = AvgTimeGraph(histTransp1_t1seed, nSteps, transpMin, transpMax, "Transp1_t1seed", "", "Seed 1 Transparency", "Average t_{1} Seed");
    TF1 *fit2_tseed = AvgTimeGraph(histTransp2_t2seed, nSteps, transpMin, transpMax, "Transp2_t2seed", "", "Seed 2 Transparency", "Average t_{2} Seed");
    AvgTimeGraph(histTransp1_t1rawseed, nSteps, transpMin, transpMax, "Transp1_t1rawseed", "", "Seed 1 Transparency", "Average t_{1} Raw Seed");
    AvgTimeGraph(histTransp2_t2rawseed, nSteps, transpMin, transpMax, "Transp2_t2rawseed", "", "Seed 2 Transparency", "Average t_{2} Raw Seed");
    AvgTimeGraph(histTransp1_t1calibseed, nSteps, transpMin, transpMax, "Transp1_t1calibseed", "", "Seed 1 Transparency", "Average t_{1} Calib Seed");
    AvgTimeGraph(histTransp2_t2calibseed, nSteps, transpMin, transpMax, "Transp2_t2calibseed", "", "Seed 2 Transparency", "Average t_{2} Calib Seed");
    AvgTimeGraph(histTransp1_t1calibseedsept, nSteps, transpMin, transpMax, "Transp1_t1calibseedsept", "", "Seed 1 Transparency", "Avg t_{1} Calib Seed Sept");
    AvgTimeGraph(histTransp2_t2calibseedsept, nSteps, transpMin, transpMax, "Transp2_t2calibseedsept", "", "Seed 2 Transparency", "Avg t_{2} Calib Seed Sept");

    cout<<"Generating sigma:Transparency 1:Transparency 2 Plots..."<<endl;
    hist2D( histTransp1Transp2_t, 2, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Sigma_t", "#sigma of t_{1}-t_{2}", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_tseed, 2, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Sigma_tseed", "#sigma of t_{1}-t_{2} Seed", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_trawseed, 2, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Sigma_trawseed", "#sigma of t_{1}-t_{2} Raw Seed", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_tcalibseed, 2, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Sigma_tcalibseed", "#sigma of t_{1}-t_{2} Calib Seed", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_tcalibseedsept, 2, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Sigma_tcalibseedsept", "#sigma of t_{1}-t_{2} Calib Seed Sept", "Seed 1 Transparency", "Seed 2 Transparency" );


    cout<<"Generating mean:Transparency 1:Transparency 2 Plots..."<<endl;
    hist2D( histTransp1Transp2_t, 1, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Mean_t", "Mean of t_{1}-t_{2}", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_tseed, 1, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Mean_tseed", "Mean of t_{1}-t_{2} Seed", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_trawseed, 1, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Mean_trawseed", "Mean of t_{1}-t_{2} Raw Seed", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_tcalibseed, 1, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Mean_tcalibseed", "Mean of t_{1}-t_{2} Calib Seed", "Seed 1 Transparency", "Seed 2 Transparency" );
    hist2D( histTransp1Transp2_tcalibseedsept, 1, nSteps,transpMin,transpMax, nSteps,transpMin,transpMax, "Transp1Transp2Mean_tcalibseedsept", "Mean of t_{1}-t_{2} Calib Seed Sept", "Seed 1 Transparency", "Seed 2 Transparency" );


    cout<<"Generating Events:Time 1:Transparency 1 Plots..."<<endl;
    hist2D( histTransp1t1, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp1t1", "Events", "Seed 1 Transparency", "t_{1}" );
    hist2D( histTransp1t1seed, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp1t1seed", "Events", "Seed 1 Transparency", "t_{1} Seed" );
    hist2D( histTransp1t1rawseed, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp1t1rawseed", "Events", "Seed 1 Transparency", "t_{1} Raw Seed" );
    hist2D( histTransp1t1calibseed, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp1t1calibseed", "Events", "Seed 1 Transparency", "t_{1} Calib Seed" );
    hist2D( histTransp1t1calibseedsept, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp1t1calibseedsept", "Events", "Seed 1 Transparency", "t_{1} Calib Seed Sept" );


    cout<<"Generating Events:Time 2:Transparency 2 Plots..."<<endl;
    hist2D( histTransp2t2, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp2t2", "Events", "Seed 2 Transparency", "t_{2}" );
    hist2D( histTransp2t2seed, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp2t2seed", "Events", "Seed 2 Transparency", "t_{2} Seed" );
    hist2D( histTransp2t2rawseed, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp2t2rawseed", "Events", "Seed 2 Transparency", "t_{2} Raw Seed" );
    hist2D( histTransp2t2calibseed, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp2t2calibseed", "Events", "Seed 2 Transparency", "t_{2} Calib Seed" );
    hist2D( histTransp2t2calibseedsept, 0, nSteps,transpMin,transpMax, nSteps,tMin,tMax, "events_Transp2t2calibseedsept", "Events", "Seed 2 Transparency", "t_{2} Calib Seed Sept" );
  }

  // THE FOLLOWING SECTION HAD IMPLEMENTED A TRANSPARENCY CORRECTION, BUT IT WAS NEGLIGIBLE SO IT WAS PHASED OUT
  /*  
  tree->SetBranchStatus("ecalElectronRechit_E",0);
  tree->SetBranchStatus("ecalElectronRechit_Eta",0);
  tree->SetBranchStatus("ecalElectronRechit_Phi",0);
  tree->SetBranchStatus("ecalElectronRechit_transpCorr",0);
  tree->SetBranchStatus("ecalElectronRechit_calibT_lagacy",0);

  tree->SetBranchStatus("t1raw_seed",0);
  tree->SetBranchStatus("t1calib_seed",0);
  tree->SetBranchStatus("t1calib_seed_sept",0);
  tree->SetBranchStatus("ele1SeedIEta",0);
  tree->SetBranchStatus("ele1SeedIPhi",0);
  tree->SetBranchStatus("ele1Eta",0);
  tree->SetBranchStatus("ele1Phi",0);

  tree->SetBranchStatus("t2raw_seed",0);
  tree->SetBranchStatus("t2calib_seed",0);
  tree->SetBranchStatus("t2calib_seed_sept",0);
  tree->SetBranchStatus("ele2SeedIEta",0);
  tree->SetBranchStatus("ele2SeedIPhi",0);
  tree->SetBranchStatus("ele2Eta",0);
  tree->SetBranchStatus("ele2Phi",0);

  cout<<"Looping through events again, subtracting avg time for transparency."<<endl;
  TH1F *histRun_tnew[nSteps];
  TH1F *histRun_tseednew[nSteps];
  TH1F *check[nSteps];
  for (int i=0; i<nSteps; i++) {
    histRun_tnew[i] = new TH1F( Form("histRun_tnew[%d]",i),";t_{1}-t_{2};Entries", 120, -3, 3);
    histRun_tseednew[i] = new TH1F( Form("histRun_tseednew[%d]",i),";t_{1}-t_{2} Seed;Entries", 120, -3, 3);
    check[i] = new TH1F( Form("check[%d]",i),"",12000, -3, 3);
  }

  for (Long64_t iEntry=0; iEntry<nentries; iEntry++) {
    if (iEntry %500000 == 0) cout << "Processing Event " << iEntry << "\n";
    tree->GetEntry(iEntry);

    if( !(ele1Pt>30 && ele2Pt>30 && mass>75 && mass<105 && ele1IsEB && ele2IsEB) ) continue; //Prelim cuts for everything


    //Time resolution vs Run
    for (int i=0; i<nSteps; i++) {
      if( !(run>=runCuts[i] && run<runCuts[i+1]) ) continue;
      float t1new, t1seednew, t2new, t2seednew;
      // Update t1
      for (int j=0; j<nSteps; j++) {
        if ( !(seed1_transpCorr>=transpCuts[j] && seed1_transpCorr<transpCuts[j+1]) ) continue;
          t1new = t1 - histTransp1_t1[j]->GetMean(); 
          t1seednew = t1_seed - histTransp1_t1seed[j]->GetMean();
          check[j]->Fill( t1seednew );
          break;
        
      }
      // Update t2
      for (int j=0; j<nSteps; j++) {
        if ( seed2_transpCorr>=transpCuts[j] && seed2_transpCorr<transpCuts[j+1] ) {
          t2new = t2 - histTransp2_t2[j]->GetMean();
          t2seednew = t2_seed - histTransp2_t2seed[j]->GetMean();
          break;
        }
      }
      // Fill histograms
      histRun_tnew[i]->Fill( t1new-t2new );
      histRun_tseednew[i]->Fill( t1seednew-t2seednew );
      break; // Don't bother to check other cases in the for loop
    }
  }
    
  SigMeanTGraph(file, histRun_tnew, "histRun_tnew", runCuts, stepSize, nSteps, "Corrected", "Run Number", "run_tnew");
  SigMeanTGraph(file, histRun_tseednew, "histRun_tseednew", runCuts, stepSize, nSteps, "Seed Corrected", "Run Number", "run_tseednew");
  AvgTimeGraph(check, nSteps, transpMin, transpMax, "Transp1_check", "", "Seed 1 Transparency", "Average t_{1} Seed Corrected");
*/


  // For reference, here are all the variables that have been activated:
  /*
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("run",1);
  tree->SetBranchStatus("nPV",1);
  tree->SetBranchStatus("eventTime",1);
  tree->SetBranchStatus("mass",1);
  tree->SetBranchStatus("ecalElectronRechit_E",1);
  tree->SetBranchStatus("ecalElectronRechit_Eta",1);
  tree->SetBranchStatus("ecalElectronRechit_Phi",1);
  tree->SetBranchStatus("ecalElectronRechit_transpCorr",1);
  tree->SetBranchStatus("ecalElectronRechit_calibT_lagacy",1);

  tree->SetBranchStatus("t1",1);
  tree->SetBranchStatus("t1_seed",1);
  tree->SetBranchStatus("t1raw_seed",1);
  tree->SetBranchStatus("t1calib_seed",1);
  tree->SetBranchStatus("t1calib_seed_sept",1);
  tree->SetBranchStatus("ele1SeedIEta",1);
  tree->SetBranchStatus("ele1SeedIPhi",1);
  tree->SetBranchStatus("ele1Eta",1);
  tree->SetBranchStatus("ele1Phi",1);
  tree->SetBranchStatus("ele1IsEB",1);
  tree->SetBranchStatus("ele1Pt",1);
  tree->SetBranchStatus("seed1_transpCorr",1);

  tree->SetBranchStatus("t2",1);
  tree->SetBranchStatus("t2_seed",1);
  tree->SetBranchStatus("t2raw_seed",1);
  tree->SetBranchStatus("t2calib_seed",1);
  tree->SetBranchStatus("t2calib_seed_sept",1);
  tree->SetBranchStatus("ele2SeedIEta",1);
  tree->SetBranchStatus("ele2SeedIPhi",1);
  tree->SetBranchStatus("ele2Eta",1);
  tree->SetBranchStatus("ele2Phi",1);
  tree->SetBranchStatus("ele2IsEB",1);
  tree->SetBranchStatus("ele2Pt",1);
  tree->SetBranchStatus("seed2_transpCorr",1);
  */
}
