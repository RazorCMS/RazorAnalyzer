RazorAnalyzer
=============

Class for analyzing the 2015 razor ntuples

Setup
-------------
From t3-higgs:

    git clone https://github.com/RazorCMS/RazorAnalyzer.git
    cd RazorAnalyzer
    source compile.sh
  
In the future, the simple compile script should be replaced with a Makefile.  

Defining a new analysis
-------------
1) Copy src/DummyAnalysis.cc and modify it to define your analyzer's behavior.  Be sure to change the name of the DummyAnalysis() function.

2) Add your analysis under "LIST OF ANALYSES" in include/RazorAnalyzer.h.

3) In src/RazorRunner.cc, under "EXECUTE YOUR ANALYSIS" add the option to execute your analysis code.

Running
------------
After compiling, 

    ./RazorRunner <list of input files> <analysis type>
  
Example: to execute a dummy analysis that does nothing,

    ./RazorRunner lists/TTJets_List_Test.txt dummy
