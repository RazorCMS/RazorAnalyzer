#!/bin/sh
g++ src/*.cc -o RazorRunner -Iinclude `root-config --cflags` `root-config --glibs`
