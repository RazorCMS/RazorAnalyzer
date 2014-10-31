include Makefile.inc

DIRS = postprocessing
SRCDIR = src
INCLUDEDIR = include

AUX = $(wildcard $(SRCDIR)/RazorAux*.cc)
ANALYSES = $(wildcard $(SRCDIR)/analyses/*.cc)

all: RazorRun
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) ); done
clean:
	@rm -f $(SRCDIR)/*.o
	@for d in $(DIRS); do (cd $$d; $(MAKE) $(MFLAGS) clean ); done

$(SRCDIR)/RazorEvents.o: $(SRCDIR)/RazorEvents.C
	$(CXX) $(SRCDIR)/RazorEvents.C $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

$(SRCDIR)/RazorAnalyzer.o: $(SRCDIR)/RazorAnalyzer.cc
	$(CXX) $(SRCDIR)/RazorAnalyzer.cc $(CXXFLAGS) -I$(INCLUDEDIR) -c $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)

RazorRun: $(SRCDIR)/RazorEvents.o $(SRCDIR)/RazorAnalyzer.o $(ANALYSES) $(AUX) $(SRCDIR)/RazorRun.cc
	$(CXX) $(SRCDIR)/RazorRun.cc $(SRCDIR)/RazorEvents.o $(ANALYSES) $(AUX) $(SRCDIR)/RazorAnalyzer.o $(CXXFLAGS) -I$(INCLUDEDIR) $(LDFLAGS) $(LIBS) -o $@ $(CXX11FLAGS)
