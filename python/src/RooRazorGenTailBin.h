//---------------------------------------------------------------------------
#ifndef ROO_RazorGenTailBin
#define ROO_RazorGenTailBin
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
class RooRazorGenTailBin : public RooAbsPdf
{
public:
   RooRazorGenTailBin() {} ;
   RooRazorGenTailBin(const char *name, const char *title,
		    RooAbsReal& _th1x,
		    RooAbsReal& _x0, RooAbsReal& _y0,
		    RooAbsReal& _b, RooAbsReal& _n,
		    RooAbsReal& _y1, RooAbsReal& _y2,
		    RooAbsReal& _xCut, RooAbsReal& _yCut, RooAbsReal& _zCut);
     //TH3* _Hnominal);
   RooRazorGenTailBin(const RooRazorGenTailBin& other,
      const char* name = 0);
   void setTH3Binning(TH3* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooRazorGenTailBin(*this,newname); }
   inline virtual ~RooRazorGenTailBin() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:
   
   RooRealProxy th1x;        // dependent variable
   RooRealProxy BX;       // X offset
   RooRealProxy BY;       // Y offset
   RooRealProxy BXX;        // shape parameter
   RooRealProxy BXY;        // shape parameter
   RooRealProxy BYY;        // shape parameter
   RooRealProxy ALPHA;        // shape parameter
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
  ClassDef(RooRazorGenTailBin,1) // RooRazorGenTailBin function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class RazorFunctionGenTail: public ROOT::Math::IParametricFunctionMultiDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(const double* x,const double* p) const
   {     
     double myexp = pow(p[0]*x[0]/1.E+3 + p[1]*x[1] + p[2]*x[0]*x[0]/1.E+6 + p[3]*x[0]*x[1]/1.E+3 + p[4]*x[1]*x[1], p[5]);
     return exp(-myexp);
   }
   
   double DoEval(const double* x) const
   {
     double myexp = pow(pars[0]*x[0]/1.E+3 + pars[1]*x[1] + pars[2]*x[0]*x[0]/1.E+6 + pars[3]*x[0]*x[1]/1.E+3 + pars[4]*x[1]*x[1], pars[5]);
     return exp(-myexp);
   }
 
   unsigned int NDim() const{
      return 2;
   }
 
   ROOT::Math::IParametricFunctionMultiDim* Clone() const
   {
      return new RazorFunctionGenTail();
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
