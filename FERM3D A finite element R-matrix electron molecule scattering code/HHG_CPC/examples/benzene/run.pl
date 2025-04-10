#!/usr/bin/perl
# This script runs the example code for benzene
# 1) run the Gaussian script "benzene.cube.inp" 
# or retrieve the potential files at
# http://fermion.colorado.edu/~tonzani/Software/rmatrix/Software.html

#2) Run gridding utility
$BIN_DIR=$PWD . "/../../bin";
system "$BIN_DIR/CM.x";
system "mv fort.16 input_grid.dat";

#3) run main code
system "$BIN_DIR/FERM3D.x";
