include Makefile.inc

DIRS = postprocessing
SRCDIR = src
INCLUDEDIR = include
INCLUDELIST= SimpleTable.h Linkdef.h

AUX = $(wildcard $(SRCDIR)/RazorAux*.cc)
ANALYSES = $(wildcard $(SRCDIR)/analyses/*.cc)

all: RazorRun NormalizeNtuple
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) ); done
clean:
	@rm -f $(SRCDIR)/*.o
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

$(INCLUDEDIR)/rootdict.cc:
	$(ROOTSYS)/bin/rootcint -f $@ -c $(CINTINCLUDES) -I$(INCLUDEDIR) $(INCLUDELIST)

$(SRCDIR)/RazorEvents.o: $(SRCDIR)/RazorEvents.C
	$(CXX) $(SRCDIR)/RazorEvents.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

RazorRun: $(SRCDIR)/RazorEvents.o $(SRCDIR)/RazorAnalyzer.o $(ANALYSES) $(AUX) $(SRCDIR)/RazorRun.cc
	$(CXX) $(SRCDIR)/RazorRun.cc $(SRCDIR)/RazorEvents.o $(ANALYSES) $(AUX) $(SRCDIR)/RazorAnalyzer.o $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

NormalizeNtuple: $(SRCDIR)/SimpleTable.o $(SRCDIR)/NormalizeNtuple.cc $(INCLUDEDIR)/rootdict.o
	$(CXX) $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(SRCDIR)/SimpleTable.o: $(SRCDIR)/SimpleTable.cc 
	$(CXX) -c $^ $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

