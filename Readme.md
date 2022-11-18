# RDF and Energy Weighting Code

This is a code to undertake calculations of derivatives of the radial distribution function, using a fluctuation theory approach.

## Dependencies 

This package uses aspects of python from the following packages, all of which must be installed prior to first use of this code.

1) Numpy

2) MDAnalysis

3) Matplotlib

## Installation Instructions

Installation should be simple, just modify the path to where you place module files on your cluster in the Makefile (the MODLOC variable).

Once you do this, then just type make!

## Contact

Copyright November 2022, Boston University

For questions, contact Zeke Piskulich (piskulichz@gmail.com)

# setup_rdf.py

This code writes input files for [LAMMPS] that calculates the energies of various groups, decided by a distance cutoff.

The essential idea of this calculation is to take the energy, H, and split it into interactions such that you get components from solute-solute, solute-close, solute-far, close-close, close-far, and far-far. Things work well when you choose a cutoff that makes close-far, far-far, and solute-far uncorrelated from the RDF. 
