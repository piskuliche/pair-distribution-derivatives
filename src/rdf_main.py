#!/usr/bin/env python

import numpy as np
import argparse
import pickle
import math
import MDAnalysis as mda

from rdf_class import RDFS
from MDAnalysis.analysis.leaflet import LeafletFinder



def Main(Iargs):
    def LIPIDS_Main(Iargs):
        u = mda.Universe(Iargs.data, Iargs.trj)
        all_rdfs = RDFS(dr = Iargs.dr, rmax = Iargs.rmax, dim = Iargs.dim)

        start, stop = Iargs.seg*Iargs.fcount, (Iargs.seg+1)*Iargs.fcount
        count = 0
        for ts in u.trajectory[start:stop:Iargs.skip]:
            leafs = LeafletFinder(u, Iargs.leafsel, pbc=True)
            leaf1, leaf2 = leafs.groups(0), leafs.groups(1)
            if count == 0: 
                all_rdfs._calc_rdf(leaf1.positions, ts, leaf1.resids, init=True)
                all_rdfs._calc_rdf(leaf2.positions, ts, leaf2.resids, init=True)
            else:
                all_rdfs._calc_rdf(leaf1.positions, ts, leaf1.resids, init=False)
                all_rdfs._calc_rdf(leaf2.positions, ts, leaf2.resids, init=False)
            count += 1
        pickle.dump(all_rdfs,open("all_rdfs.pckl",'wb'))
        return
    def REG_Main(Iargs):
        u = mda.Universe(Iargs.data, Iargs.trj)
        sel1 = u.select_atoms("Iargs.sel1")
        all_rdfs = RDFS(dr = Iargs.dr, rmax = Iargs.rmax, dim = Iargs.dim)

        start, stop = Iargs.seg*Iargs.fcount, (Iargs.seg+1)*Iargs.fcount
        count = 0
        for ts in u.trajectory[start:stop:Iargs.skip]:
            if count == 0: 
                all_rdfs._calc_rdf(sel1.positions, ts, sel1.resids, init=True)
            else:
                all_rdfs._calc_rdf(sel1.positions, ts, sel1.resids, init=False)
            count += 1
        pickle.dump(all_rdfs,open("all_rdfs.pckl",'wb'))
        return




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate Radial Distribution Function')

    parser.add_argument('-fcount',  default=10000,          type=int,   help = 'Number of fps')
    parser.add_argument('-dim',     default=2,              type=int,   help = 'Number of dimensions')
    parser.add_argument('-dr',      default=0.1,            type=float, help = 'Bin thickness')
    parser.add_argument('-rmax',    default=12,             type=float, help = 'RDF cutoff')
    parser.add_argument('-nmol',    default=288,            type=int,   help = 'Number of molecules')
    parser.add_argument('-data',    default="equil.data",   type=str,   help = 'Data file')
    parser.add_argument('-trj',     default="dump.lammpsdump",type=str, help = 'Dump file')
    parser.add_argument('-leafsel', default="type 4",       type=str,   help = 'Leaflet Selection Text')
    parser.add_argument('-seg',     default=0,              type=int,   help = 'Which segment to do?')
    parser.add_argument('-skip',    default=10,             type=int,   help = 'Frames to skip')
    parser.add_argument('-software',default='LAMMPS',       type=str,   help = 'MD Software Package')
    parser.add_argument('-lipids',  default=False,          type=bool,  help = 'True or False, incorporates leaflet info')
    
    Iargs = parser.parse_args()

    Main(Iargs)

