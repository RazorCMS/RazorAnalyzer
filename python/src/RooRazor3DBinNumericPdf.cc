//---------------------------------------------------------------------------
#include "RooFit.h"

#include "Riostream.h"
#include <TMath.h>
#include <cassert>
#include <cmath>
#include <math.h>

#include "RooRazor3DBinNumericPdf.h"
#include "RooRealVar.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/GSLIntegrator.h"
#include "Math/IntegratorMultiDim.h"

using namespace std;

ClassImp(RooRazor3DBinNumericPdf)
//---------------------------------------------------------------------------
RooRazor3DBinNumericPdf::RooRazor3DBinNumericPdf(
        const char *name, const char *title, 
        RooAbsReal& _th1x, RooAbsReal &_p0, 
        RooAbsReal &_p1, RooAbsReal &_p2, 
        RooAbsReal &_p3, RooAbsReal& _xCut, 
        RooAbsReal& _yCut, RooAbsReal& _zCut
        ) : RooAbsPdf(name, title), 
  th1x("th1x", "th1x Observable", this, _th1x),
  p0("p0", "parameter 0", this, _p0),
  p1("p1", "parameter 1", this, _p1),
  p2("p2", "parameter 2", this, _p2),
  p3("p3", "parameter 3", this, _p3),
  xCut("xCut", "X Cut parameter",this, _xCut),
  yCut("yCut", "Y Cut parameter",this, _yCut),
  zCut("zCut", "Z Cut parameter",this, _zCut),
  xBins(0),
  yBins(0),
  zBins(0),
  xMax(0),
  yMax(0),
  zMax(0),
  xMin(0),
  yMin(0),
  zMin(0),
  relTol(1E-12),
  absTol(1E-12)
{
  memset(&xArray, 0, sizeof(xArray));
  memset(&yArray, 0, sizeof(yArray));
  memset(&zArray, 0, sizeof(zArray));
}
//---------------------------------------------------------------------------
RooRazor3DBinNumericPdf::RooRazor3DBinNumericPdf(const RooRazor3DBinNumericPdf& other, const char* name) :
   RooAbsPdf(other, name), 
   th1x("th1x", this, other.th1x),   
   p0("p0", this, other.p0),
   p1("p1", this, other.p1),
   p2("p2", this, other.p2),
   p3("p3", this, other.p3),
   xCut("xCut", this, other.xCut),
   yCut("yCut", this, other.yCut),
   zCut("zCut", this, other.zCut),
   xBins(other.xBins),
   yBins(other.yBins),
   zBins(other.zBins),
   xMax(other.xMax),
   yMax(other.yMax),
   zMax(other.zMax),
   xMin(other.xMin),
   yMin(other.yMin),
   zMin(other.zMin),
   relTol(other.relTol),
   absTol(other.absTol)
{
  //memset(&xArray, 0, sizeof(xArray));
  //memset(&yArray, 0, sizeof(yArray));
  //memset(&zArray, 0, sizeof(zArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] = other.xArray[i];
  }
  for (Int_t j=0; j<yBins+1; j++){
    yArray[j] =  other.yArray[j];
  }
  for (Int_t k=0; k<zBins+1; k++){
    zArray[k] =  other.zArray[k];
  }
}
//---------------------------------------------------------------------------
void RooRazor3DBinNumericPdf::setTH3Binning(TH3* _Hnominal){
  xBins = _Hnominal->GetXaxis()->GetNbins();
  yBins = _Hnominal->GetYaxis()->GetNbins();
  zBins = _Hnominal->GetZaxis()->GetNbins();
  xMin = _Hnominal->GetXaxis()->GetBinLowEdge(1);
  yMin = _Hnominal->GetYaxis()->GetBinLowEdge(1);
  zMin = _Hnominal->GetZaxis()->GetBinLowEdge(1);
  xMax = _Hnominal->GetXaxis()->GetBinUpEdge(xBins);
  yMax = _Hnominal->GetYaxis()->GetBinUpEdge(yBins);
  zMax = _Hnominal->GetZaxis()->GetBinUpEdge(zBins);
  memset(&xArray, 0, sizeof(xArray));
  memset(&yArray, 0, sizeof(yArray));
  memset(&zArray, 0, sizeof(zArray));
  for (Int_t i=0; i<xBins+1; i++){
    xArray[i] =  _Hnominal->GetXaxis()->GetBinLowEdge(i+1);
  }
  for (Int_t j=0; j<yBins+1; j++){
    yArray[j] =  _Hnominal->GetYaxis()->GetBinLowEdge(j+1);
  }
  for (Int_t k=0; k<zBins+1; k++){
    zArray[k] =  _Hnominal->GetZaxis()->GetBinLowEdge(k+1);
  }
}
//---------------------------------------------------------------------------
void RooRazor3DBinNumericPdf::setRelTol(double _relTol){
  relTol = _relTol;
}
//---------------------------------------------------------------------------
void RooRazor3DBinNumericPdf::setAbsTol(double _absTol){
  absTol = _absTol;
}
//---------------------------------------------------------------------------
Double_t RooRazor3DBinNumericPdf::evaluate() const
{
  Double_t integral = 0.0;

  Int_t nBins = xBins*yBins*zBins;

  Int_t iBin = (Int_t) th1x;
  if(iBin < 0 || iBin >= nBins) {
    return 0.0;
  }

  
  Int_t zBin = iBin % zBins;
  Int_t yBin = ( (iBin - zBin)/(zBins) ) % (yBins);
  Int_t xBin =  (iBin - zBin - yBin*zBins ) / (zBins*yBins);

  Double_t zLow = zArray[zBin];
  Double_t zHigh = zArray[zBin+1];
  if (zCut >= zLow and zCut < zHigh){
    Double_t xLow = xArray[xBin];
    Double_t xHigh = xArray[xBin+1];
    Double_t yLow = yArray[yBin];
    Double_t yHigh = yArray[yBin+1];
    
    // define the function to be integrated numerically
    RazorFunctionNoPrefactor func;
    double params[6];    
    params[0] = p0;
    params[1] = p1;
    params[2] = p2;
    params[3] = p3;
    params[4] = xMin;
    params[5] = yMin;
    func.SetParameters(params);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(func);
    
    double a[2] = {0,0}; // (x, y) low point
    double b[2] = {0,0}; // (x, y) high point

    if(xHigh <= xCut && yHigh <= yCut) {
      return 0.0;
    }
    else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {
      a[0] = xCut;  a[1] = yLow;
      b[0] = xHigh; b[1] = yHigh;
      integral = ig.Integral(a,b);
    }
    else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {      
      a[0] = xLow;  a[1] = yCut;
      b[0] = xHigh; b[1] = yHigh;
      integral = ig.Integral(a,b);
    }
    else if(xLow < xCut && xHigh > xCut && yLow < yCut && yHigh > yCut) {            
      a[0] = xLow;  a[1] = yLow;
      b[0] = xHigh; b[1] = yHigh;
      integral = ig.Integral(a,b);     
      a[0] = xLow;  a[1] = yLow;
      b[0] = xCut; b[1] = yCut;
      integral -= ig.Integral(a,b);
    }
    else {
      a[0] = xLow;  a[1] = yLow;
      b[0] = xHigh; b[1] = yHigh;
      integral = ig.Integral(a,b);
    }
      
  }

  return integral;

}

// //---------------------------------------------------------------------------
Int_t RooRazor3DBinNumericPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const{
  if (matchArgs(allVars, analVars, th1x)) return 1;
  return 0;
}

// //---------------------------------------------------------------------------
Double_t RooRazor3DBinNumericPdf::analyticalIntegral(Int_t code, const char* rangeName) const{

   Double_t th1xMin = th1x.min(rangeName); Double_t th1xMax = th1x.max(rangeName);
   Int_t iBinMin = (Int_t) th1xMin; Int_t iBinMax = (Int_t) th1xMax;

   Double_t integral = 0.0;
      
   Int_t nBins =  xBins*yBins*zBins;

    // define the function to be integrated numerically
    RazorFunctionNoPrefactor func;
    double params[6];    
    params[0] = p0;
    params[1] = p1;
    params[2] = p2;
    params[3] = p3;
    params[4] = xMin;
    params[5] = yMin;
    func.SetParameters(params);
    ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kADAPTIVE);
    ig.SetFunction(func);
    
    double a[2] = {0,0};
    double b[2] = {0,0};

   if(code==1) { 
     for (Int_t iBin=iBinMin; iBin<iBinMax; iBin++){
       Int_t zBin = iBin % zBins;
       Int_t yBin = ( (iBin - zBin)/(zBins) ) % (yBins);
       Int_t xBin =  (iBin - zBin - yBin*zBins ) / (zBins*yBins);
 
       Double_t zLow = zArray[zBin];
       Double_t zHigh = zArray[zBin+1];
       
       
       if(iBin < 0 || iBin >= nBins) {
	 integral += 0.0;
       }
       else{
	 if (zCut >= zLow and zCut < zHigh){
	   Double_t xLow = xArray[xBin];
	   Double_t xHigh = xArray[xBin+1];
	   Double_t yLow = yArray[yBin];
	   Double_t yHigh = yArray[yBin+1];
	   

	   if(xHigh <= xCut && yHigh <= yCut) {
	     integral += 0.0;
	   }
	   else if(xLow < xCut && xHigh > xCut && yHigh <= yCut) {
	     a[0] = xCut;  a[1] = yLow;
	     b[0] = xHigh; b[1] = yHigh;
	     integral += ig.Integral(a,b);
	   }
	   else if(yLow < yCut && yHigh > yCut && xHigh <= xCut) {      
	     a[0] = xLow;  a[1] = yCut;
	     b[0] = xHigh; b[1] = yHigh;
	     integral += ig.Integral(a,b);
	   }
	   else if(xLow < xCut && xHigh > xCut && yLow < yCut && yHigh > yCut) {            
	     a[0] = xLow;  a[1] = yLow;
	     b[0] = xHigh; b[1] = yHigh;
	     integral += ig.Integral(a,b);      
	     a[0] = xLow;  a[1] = yLow;
	     b[0] = xCut; b[1] = yCut;
	     integral -= ig.Integral(a,b);
	   }
	   else {
	     a[0] = xLow;  a[1] = yLow;
	     b[0] = xHigh; b[1] = yHigh;
	     integral += ig.Integral(a,b);
	   }
	 }
       }
     }
   } else {
     cout << "WARNING IN RooRazor3DBinNumericPdf: integration code is not correct" << endl;
     cout << "                           what are you integrating on?" << endl;
     return 1.0;
   }

   return integral;
}
// //---------------------------------------------------------------------------

