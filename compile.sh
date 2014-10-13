#!/bin/sh
g++ -c src/RazorEvents.C -o src/RazorEvents.o -Iinclude `root-config --cflags` `root-config --glibs`
g++ -c src/RazorAnalyzer.cc -o src/RazorAnalyzer.o -Iinclude `root-config --cflags` `root-config --glibs`
g++ src/RazorRunner.cc src/analyses/*.cc src/*.o -o RazorRunner -Iinclude `root-config --cflags` `root-config --glibs`
