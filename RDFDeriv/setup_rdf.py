#!/usr/bin/env python
import glob, shutil
import numpy as np
import MDAnalysis as mda
import sys


"""

Copyright November 2022 - Boston University

This is a code that sets up an energy calculation for different groups, split up into "solute", "close", and "far" for each molecule in the system. It right now works for primarily lammps, but my goals is to extend this eventually to gromacs as well.

"""

def LAMMPS_Main(datafile, trjfile, sysfile, nconfigs=5000, skip=1, framesep=1000,
                 startframe=1000000, dim=3, cut=12.0, subdir="rdf_test"):
    """Calls the code for writing LAMMPS Files

    This function calls the code for writing LAMMPS files for calculating the energies.

    Logfiles are written to subdir_setup.log

    Args:
        datafile (str): Filename of lammps data file
        trjfile (str): Filename of trajectory file
        sysfile (str): Filename of system input file with run commands removed 
        nconfigs (int): Number of configurations [default=5000]
        skip (int): Number of configs to skip [default=1]
        framesep (int): Number of timesteps between frames [default=1000]
        startframe (int): Starting frame number [default-1000000]
        dim (int): Number of dimensions [default=2]
        cut (float): Spherical/cylindrical cutoff values [default=12.0]
        subdir (str): Name of the subdirectory to generate the files. [default=rdf_test]
        

    """
    u = mda.Universe(datafile, trjfile)
    atoms = u.select_atoms("all")
    nmols=0
    logger = open(subdir+"_setup.log",'w')
    for ts in u.trajectory[0:nconfigs:skip]:
        current_frame = ts.frame*framesep + startframe
        logger.write("Working on frame: %d\n" % current_frame)
        L = ts.dimensions[:3]
        comr = atoms.center_of_mass(compound='residues')
        nmols = np.shape(comr)[0]
        if ts.frame == 0:
            Prep_Dir(nmols, subdir)
        dr = Get_Distances(comr,L,dim)
        Write_Groups(dr, current_frame, subdir, cut)
    Write_System(nmols, sysfile=sysfile, subdir=subdir)
    Write_Task(nmols = nmols, queue_engine="TORQUE", nmol_per_task=100, hours=2, procs=4, subdir=subdir)
    logger.close()

def GMX_Main(datafile, trjfile, nconfigs=5000, skip=1, framesep=1000, startframe=1000000, dim=3):
    """Calls the code for writing Gromacs Files

    This function calls the code for writing GROMACS files for calculating the energies.

    Args:
        Iargs (ArgParse Object): These are the system input arguments used by argparse.
    Todo:
        * Implement these features for gromacs, may potentially require writing some new functions.
    """
    exit("Error: Gromacs support not yet added")

def Get_Distances(r,L,dim=3):
    """Calculates the self distance array between a vector and its other elements.

    This uses numpy broadcasting to calculate ALL the pairwise distances and outputs them in a matrix of distances.

    This is described here: 
    https://stackoverflow.com/questions/60039982/numpy-python-vectorize-distance-function-to-calculate-pairwise-distance-of-2-ma/60040269#60040269
    
    Args:
        r (array_like): Array of positions shape(natoms,dim)
        L (array_like): Box dimension array shape(dim,)
        dim (int): Box dimensionality [default=3]

    Returns:
        array_like: returns pairwise distances

    """
    vecdr = r[:,np.newaxis,:dim]-r[np.newaxis,:,:dim]
    vecdr = vecdr - np.multiply(L[:dim],np.round(np.divide(vecdr,L[:dim])))
    dr = np.linalg.norm(vecdr,axis=-1)
    return dr

def Sort_Distances(dr):
    """ This function takes an array (dr) and sorts it by the last axis.

    It returns:
    1) the idx map of the sort
    2) the sorted dr array (sdr)

    Args: 
        dr (array_like): Pairwise distances array
    Returns:
        int: index map of the sorted array
        array_like: Sorted array

    """
    natoms = np.shape(dr)[0]
    idx = np.argsort(dr,axis=-1,kind='quicksort')
    sdr=np.take_along_axis(dr, idx, axis=-1)
    return idx,sdr

def Prep_Dir(nmols, subdir):
    """Write files and directories to compute energies for LAMMPS
    
    This function prepares files to calculate the energies for each group (solute, close, and far) to 
    get the energies. This makes any directories that do not already exist.

    This writes files in includ_dir, and calc_dir, which are subdirs of Iargs.subdir.

    For lammps, we use pe/tally because it is faster than compute group/group. 

    Args:
        nmols (int): Number of molecules
        subdir (str): Name of subdirectory for calculation.
    
    """
    import os
    if not os.path.exists("./%s" % subdir):
        os.makedirs("./%s" % subdir)
    if not os.path.exists("./%s/include_dir" % subdir):
        os.makedirs("./%s/include_dir" % subdir)
    if not os.path.exists("./%s/ener_dir" % subdir):
        os.makedirs("./%s/ener_dir" % subdir)
    if not os.path.exists("./%s/calc_dir" % subdir):
        os.makedirs("./%s/calc_dir" % subdir)
    if not os.path.exists("./%s/calc_dir/logs" % subdir):
        os.makedirs("./%s/calc_dir/logs" % subdir)
    for basefile in glob.glob("./base/*"):
        shutil.copy(basefile,"./%s/" % subdir)

    fi = open("./%s/include_dir/include.computes" % subdir,'w')
    # Define the initial groups to be 0
    for grp in ["lt","solu","close","far"]:
        fi.write("group %s molecule 0\n"%grp)

    # Pairwise Components
    fi.write("compute sspair solu group/group solu kspace yes pair yes\n")
    fi.write("compute ccpair close group/group close kspace yes pair yes\n")
    fi.write("compute ffpair far group/group far kspace yes pair yes\n")
    fi.write("compute scpair solu group/group close kspace yes pair yes\n")
    fi.write("compute sfpair solu group/group far kspace yes pair yes\n")
    fi.write("compute cfpair close group/group far kspace yes pair yes\n")


    # Per Atom components
    fi.write("compute sintrapa solu pe/atom bond angle\n")
    fi.write("compute cintrapa close pe/atom bond angle\n")
    fi.write("compute fintrapa far pe/atom bond angle\n")
    # Sum those components
    fi.write("compute sintra solu reduce sum c_sintrapa\n")
    fi.write("compute cintra close reduce sum c_cintrapa\n")
    fi.write("compute fintra far reduce sum c_fintrapa\n")
    fi.close()
        
    for i in range(nmols):
        incf = open("./%s/include_dir/include.groups-%d"%(subdir,i),'w')
        # Final Files
        incf.write("fix wpair all ave/time 1000 1 1000 c_sspair c_ccpair c_ffpair c_scpair c_sfpair c_cfpair file ../ener_dir/pair.compute%d\n"%i)
        incf.write("fix wintr all ave/time 1000 1 1000 c_sintra c_cintra c_fintra file ../ener_dir/intra.compute%d\n"%i)
        incf.close()

def Write_Groups(dr, frame, subdir, cut=12):
    """Writes occupancies in each group to a lammps file
    
    This writes the group occupancies out to a different file for each molcule. This is based on a distance cutoff,
    which is stored in Iargs.cut

    Args:
        dr (array_like): Array of pairwise distances
        frame (int): Timestep associated with the current frame
        subdir (str): Subdirectory where rdf calculation will take place
        cut (float): Spherical (or cylindrical) cutoff distance for defining close group

    """

    nmols=np.shape(dr)[0]
    mols = np.arange(0,nmols)
    for i in range(nmols):
        fi = open("./%s/include_dir/include.groups-%d"%(subdir, i),'a')
        close=np.where(dr[i]<cut)[-1]
        for grp in ["lt","solu","close","far"]:
            fi.write("group %s clear\n"%grp)
        fi.write("\n")
        fi.write("group close molecule ")
        for j in close:
            if j != i: fi.write("%d "% (j+1)) # Ignore i from close
        fi.write("\n")
        fi.write("group solu molecule %d\n"%(i+1))
        fi.write("group far subtract all close solu\n")

        fi.write("\n")
        fi.write('print "beginning rerun %d"\n'%frame)
        fi.write("rerun ../../dump_dir/dump.lammpsdump-%d first %d last %d dump x y z\n" % (frame,frame,frame))
        fi.write('print "finishing rerun"\n')
        fi.close()

def Write_System(nmols, sysfile="system.in", subdir="rdf_test"):
    """This writes the includes into the simulation input file.

    This takes a base system input file (Iargs.sys) and writes them into a separate system input file
    for each molecule, so that lammps can easily calculate all of them.

    Args:
        nmols (int): Number of molecules
        sysfile (str): System input file with no run commands inside of it.
        subdir (str): Subdirectory where rdf calculation will take place

    """
    print("Writing System Files")
    lines=None
    with open("./%s/%s" %(subdir, sysfile),'r') as f:
        lines = f.readlines()
    for i in range(nmols):
        with open("./%s/calc_dir/%s-%d"%(subdir,sysfile,i),'w') as g:
            for line in lines:
                g.write(line)
            g.write("log logs/log.%d.out\n"%i)
            g.write("include ../include_dir/include.computes\n")
            g.write("include ../include_dir/include.groups-%d\n"%i)
    
def Write_Task(nmols = 343, queue_engine="TORQUE", nmol_per_task=1, hours=2, procs=4, subdir="rdf_test"):
    """This is a simple function to write a submit script

    Writes a file, called submit_task.sh

    Args:
        nmols (int): Number of molecules
        queue_engine (str): What queue engine should it write for?
            Options are TORQUE, or SLURM
        nmol_per_task (int): Number of molecules per task
        hours (int): Number of calculation hours
        procs (int): Number of processors (must be appropriate for system)
        subdir (str): Subdirectory where rdf calculation will take place

    Todo:
        * SLURM Support

    """
    try:
        assert type(nmol_per_task) == int
    except:
        raise AssertionError("Error: nmol_per_task must be an integer")
    
    try:
        assert queue_engine in ["SLURM","TORQUE"]
    except:
        raise AssertionError("Error: queue_engine %s is not yet supported" % queue_engine)
    import numpy as np

    stop_array = np.ceil(nmols/nmol_per_task)
    
    if queue_engine == "TORQUE":
        with open("./%s/submit_task.sh" % subdir,'w') as f:
            f.write("#!/bin/bash\n")
            f.write("#\n")
            f.write("#$ -l h_rt=%d:00:00\n"%hours)
            f.write("#$ -j y\n")
            f.write("#$ -N ener-task\n")
            f.write("#$ -pe omp %d\n" % procs)
            f.write("#$ -V\n")
            if stop_array > 1: f.write("#$ -t 1-%d\n" % stop_array)
            f.write("\n")

            if stop_array > 1: 
                f.write("fid=$((SGE_TASK_ID-1))\n\n")
            else:
                f.write("fid=0")

            f.write("module load gcc/8.3.0\n")
            f.write("module load openmpi/3.1.4_gnu-8.1\n")
            f.write("module load fftw/3.3.8\n")
            f.write("module load lammps/3Mar2020\n")
            f.write("\n")

            f.write("cd calc_dir\n")
            f.write("for i in {1..%d}; do\n" % nmol_per_task)
            f.write("   molid=$((fid*%d + i - 1))\n" % nmol_per_task)
            f.write("   if [ ${molid} -gt %d ]; then\n" % nmols)
            f.write("       break\n")
            f.write("   fi\n")
            f.write("   mpirun lmp_mpi -in system.in-${molid}\n")
            f.write("done\n")
            f.write("cd -\n")

    elif queue_engine == "SLURM":
        raise ValueError("Error: SLURM Support has not been fully added")
        with open("submit_task.sh",'w') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH -l h_rt=%d:00:00\n"%hours)
            f.write("#SBATCH -j y\n")
            f.write("#SBATCH -N ener-task\n")
            f.write("#SBATCH -pe omp %d\n" % procs)
            f.write("#SBATCH -V\n")
            if stop_array > 1: f.write("#SBATCH -t 1-%d\n" % stop_array)
            f.write("\n")

            if stop_array > 1: 
                f.write("fid=$((SLURM_ARRAY_TASK_ID-1))\n\n")
            else:
                f.write("fid=0")

            f.write("module load gcc/8.3.0\n")
            f.write("module load openmpi/3.1.4_gnu-8.1\n")
            f.write("module load fftw/3.3.8\n")
            f.write("module load lammps/3Mar2020\n")
            f.write("\n")

            f.write("cd calc_dir\n")
            f.write("for i in {1..%d}; do\n" % nmol_per_task)
            f.write("   molid=$((fid*%d + i - 1))\n" % nmol_per_task)
            f.write("   if [ ${molid} -gt %d ]; then\n" % nmols)
            f.write("       break\n")
            f.write("   fi\n")
            f.write("   mpirun lmp_mpi -in system.in-${molid}\n")
            f.write("done\n")
            f.write("cd -\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This code creates a series of files for lammps simulations that quickly rerun and re-calculate energies')
    parser.add_argument('-data',    default="equil.data",       type=str,help='Name of data file [default equil.data]')
    parser.add_argument('-dump',    default="dump.lammpsdump",  type=str,help='Name of dump file [default dump.lammpsdump]')
    parser.add_argument('-nconfigs',default=10000,              type=int, help='Number of configurations [default 10000]')
    parser.add_argument('-sys',     default="system.in",        type=str, help='System input file [default system.in]')
    parser.add_argument('-cut',     default=30,                 type=float, help='Cutoff to use for defining what is close')
    parser.add_argument('-skip',    default=1,                  type=int, help='Frequency to use frames')
    parser.add_argument('-framesep',default=1000,               type=int, help='Dump Frequency (Excluding Skip)')
    parser.add_argument('-startframe',default=1000000,          type=int, help='Starting frame timestep number')
    parser.add_argument('-dim',     default=2,                  type=int, help='Dimsensionality (2 or 3)')
    parser.add_argument('-software',default='LAMMPS',           type=str, help='MD Simualtion Program')
    parser.add_argument('-subdir',  default='rdf_deriv',        type=str, help='Subdirectory for calculation')
    Iargs = parser.parse_args()

    if Iargs.software == 'LAMMPS':
        LAMMPS_Main(Iargs.data, Iargs.dump, Iargs.sys, nconfigs=Iargs.nconfigs, skip=Iargs.skip, framesep=Iargs.framesep,
                 startframe=Iargs.startframe, dim=Iargs.dim, cut=Iargs.cut, subdir=Iargs.subdir)
    elif Iargs.software == 'GMX':
        GMX_Main(Iargs.data, Iargs.dump, Iargs.sys, nconfigs=Iargs.nconfigs, skip=Iargs.skip, framesep=Iargs.framesep,
                 startframe=Iargs.startframe, dim=Iargs.dim, cut=Iargs.cut, subdir=Iargs.subdir)
    else:
        exit("Error: Sofware %s not recognized, options are LAMMPS or GMX" % Iargs.software)