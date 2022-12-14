#!/usr/bin/env python

def LIPIDS_Main(datafile, trjfile, leafsel, dr=0.1, rmax=12.0, dim=2, fcount=5000, skip=1, segid = 0):
    """Sets up the calculation of the RDF for Lipid Membranes

    This function calculates the RDF for lipid membranes, by using MDAnalysis leaflet analysis
    and treating each of the leaflets separately, before combining them in the end. 
    Essentially it exists to call the RDFS class.

    Args:
        datafile (str): Data file name
        trjfile (str): Trajectory file name
        leafsel (str): Leaflet selection text
        dr (float): RDF bin width [default 0.1]
        rmax (float): RDF cutoff [default 12.0]
        dim (int): Dimensionality [default 2]
        fcount (int): Number of configurations per segment [default 5000]
        skip (int): Take every skipth frame [default 1]
        segid (int): Segment number [default 0]
    
    Raises:
        ValueError: Don't ask for it to use more than 50,000 configs.

    """
    u = mda.Universe(datafile, trjfile)
    all_rdfs = RDFS(dr = dr, rmax = rmax, dim = dim)

    start, stop = segid*fcount, (segid+1)*fcount
    if ((stop-start)/skip) > 50000:
        raise ValueError("Error: Too many configurations. Use the seg option to subdivide this calculation.")

    count = 0
    for ts in u.trajectory[start:stop:skip]:
        leafs = LeafletFinder(u, leafsel, pbc=True)
        leaf1, leaf2 = leafs.groups(0), leafs.groups(1)
        L = ts.dimensions[:dim]
        if count == 0: 
            all_rdfs.calc_rdf(leaf1.positions, L, leaf1.resids, init=True)
            all_rdfs.calc_rdf(leaf2.positions, L, leaf2.resids, init=True)
        else:
            all_rdfs.calc_rdf(leaf1.positions, L, leaf1.resids, init=False)
            all_rdfs.calc_rdf(leaf2.positions, L, leaf2.resids, init=False)
        count += 1
    pickle.dump(all_rdfs,open("all_rdfs-%d.pckl"%segid,'wb'))

def REG_Main(datafile, trjfile, molsel1, dr=0.1, rmax=12.0, dim=2, fcount=5000, skip=1, segid = 0):
    """Sets up the calculation of the RDF

    This function calculates the RDF for regular liquid systems. Essentially it exists to call the RDFS class.

    Args:
        datafile (str): Data file name
        trjfile (str): Trajectory file name
        leafsel (str): Leaflet selection text
        dr (float): RDF bin width [default 0.1]
        rmax (float): RDF cutoff [default 12.0]
        dim (int): Dimensionality [default 2]
        fcount (int): Number of configurations per segment [default 5000]
        skip (int): Take every skipth frame [default 1]
        segid (int): Segment number [default 0]
    
    Raises:
        ValueError: Don't ask for it to use more than 50,000 configs.

    """

    u = mda.Universe(datafile, trjfile)
    sel1 = u.select_atoms("%s"%molsel1)
    all_rdfs = RDFS(dr = dr, rmax = rmax, dim = dim)

    start, stop = segid*fcount, (segid+1)*fcount
    if ((stop-start)/skip) > 50000:
        raise ValueError("Error: Too many configurations. Use the seg option to subdivide this calculation.")

    count = 0
    for ts in u.trajectory[start:stop:skip]:
        L = ts.dimensions[:dim]
        if count == 0: 
            all_rdfs.calc_rdf(sel1.positions, L, sel1.resids, init=True)
        else:
            all_rdfs.calc_rdf(sel1.positions, L, sel1.resids, init=False)
        count += 1
    pickle.dump(all_rdfs,open("all_rdfs-%d.pckl"%segid,'wb'))
    return





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
    parser.add_argument('-sel1',    default="type 1",       type=str,   help = 'Selection')
    parser.add_argument('--lipids', dest='lipids', action='store_true')
    parser.add_argument('--no-lipids', dest='lipids', action='store_false')
    parser.set_defaults(lipids=True)

    Iargs = parser.parse_args()

    if Iargs.lipids == True:
        LIPIDS_Main(Iargs.datafile, Iargs.trjfile, Iargs.leafsel, dr=Iargs.dr, rmax=Iargs.rmax, dim=Iargs.dim, fcount=Iargs.fcount,
                    skip=Iargs.skip, segid=Iargs.seg)
    else:
        REG_Main(Iargs.datafile, Iargs.trjfile, Iargs.sel1, dr=Iargs.dr, rmax=Iargs.rmax, dim=Iargs.dim, fcount=Iargs.fcount,
                    skip=Iargs.skip, segid=Iargs.seg)

