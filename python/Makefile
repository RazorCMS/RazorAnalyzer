include ../Makefile.inc

INCLUDEDIR= src
SHARED = libRazorRun2.so
INCLUDELIST= RooRazor2DTail_SYS.h RooRazor3DBinPdf.h MLMixedSample.h RooBTagMult.h RooSigmoid.h Linkdef.h

all: $(SHARED) 

$(INCLUDEDIR)/rootdict.cc: 
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) -I$(INCLUDEDIR) $(INCLUDELIST) 

$(INCLUDEDIR)/RooRazor2DTail_SYS.o:
	$(CXX) $(INCLUDEDIR)/RooRazor2DTail_SYS.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooRazor3DBinPdf.o:
	$(CXX) $(INCLUDEDIR)/RooRazor3DBinPdf.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooBTagMult.o:
	$(CXX) $(INCLUDEDIR)/RooBTagMult.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/RooSigmoid.o:
	$(CXX) $(INCLUDEDIR)/RooSigmoid.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(INCLUDEDIR)/MLMixedSample.o:
	$(CXX) $(INCLUDEDIR)/MLMixedSample.cc $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $@ -I$(INCLUDEDIR) -c

$(SHARED): $(INCLUDEDIR)/rootdict.o $(INCLUDEDIR)/RooRazor2DTail_SYS.o $(INCLUDEDIR)/RooRazor3DBinPdf.o $(INCLUDEDIR)/RooBTagMult.o $(INCLUDEDIR)/RooSigmoid.o $(INCLUDEDIR)/MLMixedSample.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) -o $@ $^
	mv $(SHARED) lib

clean:
	@rm -f $(INCLUDEDIR)/*.o
	@rm -f $(INCLUDEDIR)/rootdict.*
	@rm -f $(SHARED) 
	@rm -f lib/$(SHARED)
