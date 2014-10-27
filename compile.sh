#!/bin/sh
g++ -c src/RazorEvents.C -o src/RazorEvents.o -Iinclude `root-config --cflags` `root-config --glibs` -std=c++0x
g++ -c src/RazorAnalyzer.cc -o src/RazorAnalyzer.o -Iinclude `root-config --cflags` `root-config --glibs` -std=c++0x
g++ src/RazorRun.cc src/analyses/*.cc src/*.o src/RazorAux*.cc -o RazorRun -Iinclude `root-config --cflags` `root-config --glibs` -std=c++0x
