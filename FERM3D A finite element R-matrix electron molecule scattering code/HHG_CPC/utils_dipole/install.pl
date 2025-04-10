#!/usr/bin/perl
# This script compiles all the subroutines and libraries
# Put compiler as variable or make multiple scripts

# Define Compilers
$FC='ifort';
#$COMPFLG='-check all -g -O0 -C';
$COMPFLG='-O2';
$CC='gcc';

open OUT, ">INSTALL.sh";
print OUT (<<"EOS");
#/usr/bin/csh
make clean
$FC $COMPFLG -c modules_3d.f90
$FC $COMPFLG -c gau3d.f90
$FC $COMPFLG -c Hydrogen_perturbed.f90

make 
make CM.x
make clean
EOS
close OUT;

system 'chmod +x INSTALL.sh';
system './INSTALL.sh';
