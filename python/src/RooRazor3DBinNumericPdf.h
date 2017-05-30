//---------------------------------------------------------------------------
#ifndef ROO_Razor3DBinNumericPdf
#define ROO_Razor3DBinNumericPdf
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH3.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"

//---------------------------------------------------------------------------
class RooRazor3DBinNumericPdf : public RooAbsPdf
{
public:
   RooRazor3DBinNumericPdf() {} ;
   RooRazor3DBinNumericPdf(const char *name, const char *title,
		    RooAbsReal& _th1x,
		    RooAbsReal& _p0, RooAbsReal& _p1,
		    RooAbsReal& _p2, RooAbsReal& _p3,
		    RooAbsReal& _xCut, RooAbsReal& _yCut, 
                    RooAbsReal& _zCut);
   RooRazor3DBinNumericPdf(const RooRazor3DBinNumericPdf& other,
      const char* name = 0);
   void setTH3Binning(TH3* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooRazor3DBinNumericPdf(*this,newname); }
   inline virtual ~RooRazor3DBinNumericPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:
   
   RooRealProxy th1x;        // dependent variable
   RooRealProxy p0;       
   RooRealProxy p1;      
   RooRealProxy p2;    
   RooRealProxy p3;   
   RooRealProxy xCut;        // X cut constant (set pdf to 0 for X < Xcut && Y < Ycut) 
   RooRealProxy yCut;        // Y cut constant (set pdf to 0 for X < Xcut && Y < Ycut) 
   RooRealProxy zCut;        // Z cut constant (set pdf to 0 unless Zut <= Z < Zcut)
   Int_t xBins;        // X bins
   Int_t yBins;        // Y bins
   Int_t zBins;        // Z bins
   Double_t xArray[20]; // xArray[xBins+1]
   Double_t yArray[20]; // yArray[yBins+1]
   Double_t zArray[5]; // zArray[zBins+1]
   Double_t xMax;        // X max
   Double_t yMax;        // Y max
   Double_t zMax;        // Z max
   Double_t xMin;        // X min
   Double_t yMin;        // Y min
   Double_t zMin;        // Z min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
  ClassDef(RooRazor3DBinNumericPdf,1) // RooRazor3DBinNumericPdf function
    
};
//---------------------------------------------------------------------------

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class RazorFunctionNoPrefactor : public ROOT::Math::IParametricFunctionMultiDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(const double* x,const double* p) const
   {     
        double xOffset = p[0] - p[4];
        double yOffset = p[1] - p[5];
        double n = p[3];
        double myexp = p[2] * pow( 
            fabs((x[0]+xOffset) * (x[1]+yOffset)),
            1/n);
       return exp(-myexp);
   }
   
   double DoEval(const double* x) const
   {
        double xOffset = pars[0] - pars[4];
        double yOffset = pars[1] - pars[5];
        double n = pars[3];
        double myexp = pars[2] * pow( 
            fabs((x[0]+xOffset) * (x[1]+yOffset)),
            1/n);
       return exp(-myexp);
   }
 
   unsigned int NDim() const{
      return 2;
   }
 
   ROOT::Math::IParametricFunctionMultiDim* Clone() const
   {
      return new RazorFunctionNoPrefactor();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
   
   unsigned int NPar() const
   {
      return 6;
   }
};

#endif
