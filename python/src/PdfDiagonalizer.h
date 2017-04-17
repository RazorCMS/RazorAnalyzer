#ifndef ROO_PdfDiagonalizer
#define ROO_PdfDiagonalizer

class RooWorkspace;
struct RooFitResult;
struct RooAbsPdf;

#include <string>
#include <RooArgList.h>

class PdfDiagonalizer {
    public:
        PdfDiagonalizer(const char *name, RooWorkspace *w, RooFitResult &result);

        RooAbsPdf *diagonalize(RooAbsPdf &pdf) ;
        const RooArgList & originalParams() { return parameters_; }
        const RooArgList & diagonalParams() { return eigenVars_; }
    private:
        std::string name_;
        RooArgList  parameters_;
        RooArgList  eigenVars_;
        RooArgList  replacements_;
};

#endif
