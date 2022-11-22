#!/usr/bin/env python





def Main(Iargs):
    # Calls the correct rdf calculation based off of input arguments.
    def LIPIDS_Main(Iargs):
        """Sets up the calculation of the RDF for Lipid Membranes

        This function calculates the RDF for lipid membranes, by using MDAnalysis leaflet analysis
        and treating each of the leaflets separately, before combining them in the end. 
        Essentially it exists to call the RDFS class.

        Args:
            Iargs (ArgParse Object): These are the system input arguments used by argparse.

        """
        u = mda.Universe(Iargs.data, Iargs.trj)
        all_rdfs = RDFS(dr = Iargs.dr, rmax = Iargs.rmax, dim = Iargs.dim)

        start, stop = Iargs.seg*Iargs.fcount, (Iargs.seg+1)*Iargs.fcount
        if ((stop-start)/Iargs.skip) > 50000:
            exit("Error: Too many configurations. Use the seg option to subdivide this calculation.")

        count = 0
        for ts in u.trajectory[start:stop:Iargs.skip]:
            leafs = LeafletFinder(u, Iargs.leafsel, pbc=True)
            leaf1, leaf2 = leafs.groups(0), leafs.groups(1)
            L = ts.dimensions[:Iargs.dim]
            if count == 0: 
                all_rdfs.calc_rdf(leaf1.positions, L, leaf1.resids, init=True)
                all_rdfs.calc_rdf(leaf2.positions, L, leaf2.resids, init=True)
            else:
                all_rdfs.calc_rdf(leaf1.positions, L, leaf1.resids, init=False)
                all_rdfs.calc_rdf(leaf2.positions, L, leaf2.resids, init=False)
            count += 1
        pickle.dump(all_rdfs,open("all_rdfs-%d.pckl"%Iargs.seg,'wb'))

    def REG_Main(Iargs):
        """Sets up the calculation of the RDF

        This function calculates the RDF for regular liquid systems. Essentially it exists to call the RDFS class.

        Args:
            Iargs (ArgParse Object): These are the system input arguments used by argparse.
            
        """
        u = mda.Universe(Iargs.data, Iargs.trj)
        sel1 = u.select_atoms("Iargs.sel1")
        all_rdfs = RDFS(dr = Iargs.dr, rmax = Iargs.rmax, dim = Iargs.dim)

        start, stop = Iargs.seg*Iargs.fcount, (Iargs.seg+1)*Iargs.fcount
        if ((stop-start)/Iargs.skip) > 50000:
            exit("Error: Too many configurations. Use the seg option to subdivide this calculation.")

        count = 0
        for ts in u.trajectory[start:stop:Iargs.skip]:
            L = ts.dimensions[:Iargs.dim]
            if count == 0: 
                all_rdfs.calc_rdf(sel1.positions, L, sel1.resids, init=True)
            else:
                all_rdfs.calc_rdf(sel1.positions, L, sel1.resids, init=False)
            count += 1
        pickle.dump(all_rdfs,open("all_rdfs-%d.pckl"%Iargs.seg,'wb'))
        return
    if Iargs.lipids == True:
        LIPIDS_Main(Iargs)
    else:
        REG_Main(Iargs)




if __name__ == "__main__":

    import numpy as np
    import argparse
    import pickle
    import math
    import MDAnalysis as mda

    from rdf_class import RDFS
    from MDAnalysis.analysis.leaflet import LeafletFinder
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
    parser.add_argument('-skip',    default=1,              type=int,   help = 'Frames to skip')
    parser.add_argument('-software',default='LAMMPS',       type=str,   help = 'MD Software Package')
    parser.add_argument('-lipids',  default=True,          type=bool,  help = 'True or False, incorporates leaflet info')
    
    Iargs = parser.parse_args()

    Main(Iargs)

