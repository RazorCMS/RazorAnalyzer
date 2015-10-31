include Makefile.inc

DIRS = python
SRCDIR = src
INCLUDEDIR = include
INCLUDELIST= SimpleTable.h Linkdef.h

SCRATCHDIR = /tmp/

AUX = $(wildcard $(SRCDIR)/RazorAux*.cc $(SRCDIR)/HggRazorAuxPhoton.cc $(SRCDIR)/Hemisphere.cc $(SRCDIR)/Davismt2.cc)
ANALYSES = $(wildcard analyses/*.cc)
JETCORR = $(SRCDIR)/JetCorrectorParameters.cc $(SRCDIR)/SimpleJetCorrectionUncertainty.cc  $(SRCDIR)/JetCorrectionUncertainty.cc $(SRCDIR)/SimpleJetCorrector.cc $(SRCDIR)/FactorizedJetCorrector.cc $(SRCDIR)/SimpleJetResolution.cc $(SRCDIR)/BTagCalibrationStandalone.cc
JETCORROBJ = $(JETCORR:cc=o)
EXECUTABLES = RazorRun NormalizeNtuple SkimNtuple

.PHONY: clean all lxplus

all: $(addprefix $(SCRATCHDIR)/, $(EXECUTABLES))
	mv $(addprefix $(SCRATCHDIR)/, $(EXECUTABLES)) .
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) ); done
lxplus: $(EXECUTABLES)
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) ); done

clean:
	@-rm $(EXECUTABLES)
	@rm -f $(SRCDIR)/*.o
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

$(INCLUDEDIR)/rootdict.cc:
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) -I$(INCLUDEDIR) $(INCLUDELIST)

$(SRCDIR)/SimpleTable.o: $(SRCDIR)/SimpleTable.cc 
	$(CXX) -c $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(SRCDIR)/RazorEvents.o: $(SRCDIR)/RazorEvents.C $(INCLUDEDIR)/RazorEvents.h
	$(CXX) $(SRCDIR)/RazorEvents.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/RazorEvents.o $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(JETCORROBJ): %.o: %.cc
	$(CXX) -c $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS) $<

$(SCRATCHDIR)/RazorRun: $(SRCDIR)/RazorEvents.o $(SRCDIR)/RazorAnalyzer.o $(ANALYSES) $(JETCORROBJ) $(AUX) $(SRCDIR)/RazorRun.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(SCRATCHDIR)/NormalizeNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/NormalizeNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(SCRATCHDIR)/SkimNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/SkimNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

RazorRun: $(SRCDIR)/RazorEvents.o $(SRCDIR)/RazorAnalyzer.o $(ANALYSES) $(JETCORROBJ) $(AUX) $(SRCDIR)/RazorRun.cc
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

NormalizeNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/NormalizeNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

SkimNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/SkimNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

MakePlots: $(SRCDIR)/SimpleTable.o ./macros/BackgroundStudies/OverlayKinematicPlots_Selected.C $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)
