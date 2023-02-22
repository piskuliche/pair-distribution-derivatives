import glob, shutil
import numpy as np

def LAMMPS_Fort_Main(molecules, nconfigs=5000, skip=1, framesep=1000, startframe=1000000, subdir="rdf_deriv", sysfile="system.in"):
    for molid in molecules:
        print("molid %d" % molid)
        neighbors = Read_Neighbor_File(molid, nconfigs=nconfigs, skip=skip, subdir=subdir)
        print(np.shape(neighbors))
        framesep = skip*framesep
        for i in range(0,nconfigs,skip):
            frame = startframe + framesep*i 
            Write_Groups(molid, neighbors[i], frame, subdir=subdir)
    Write_System(molecules, sysfile=sysfile, subdir=subdir)

    return

def Read_Neighbor_File(molid, nconfigs=5000, skip=1,subdir="rdf_deriv"):
    neighbors = []
    with open("%s/neighbors/neighbors_%d.dat"%(subdir,molid), 'r') as f:
        lines = f.readlines()
        for line in lines:
            neighbors.append(line.strip())
    return neighbors[1::skip]

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
    if not os.path.exists("./%s/neighbors" % subdir):
        raise OSError("No neighbors directory found. Please run the fortran code first.")
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
        incf = open("./%s/include_dir/include.groups-%d"%(subdir,i+1),'w')
        # Final Files
        incf.write("fix wpair all ave/time 1000 1 1000 c_sspair c_ccpair c_ffpair c_scpair c_sfpair c_cfpair file ../ener_dir/pair.compute%d\n"%i)
        incf.write("fix wintr all ave/time 1000 1 1000 c_sintra c_cintra c_fintra file ../ener_dir/intra.compute%d\n"%i)
        incf.close()

def Write_Groups(molid, neighbors, frame, subdir):
    """Writes occupancies in each group to a lammps file
    
    This writes the group occupancies out to a different file for each molcule. This is based on a distance cutoff,
    which is stored in Iargs.cut

    Args:
        dr (array_like): Array of pairwise distances
        frame (int): Timestep associated with the current frame
        subdir (str): Subdirectory where rdf calculation will take place
        cut (float): Spherical (or cylindrical) cutoff distance for defining close group

    """
    fi = open("./%s/include_dir/include.groups-%d"%(subdir, molid),'a')
    for grp in ["lt","solu","close","far"]:
        fi.write("group %s clear\n"%grp)
    fi.write("\n")
    fi.write("group close molecule %s\n"% neighbors)
    fi.write("\n")
    fi.write("group solu molecule %d\n"%(molid))
    fi.write("group far subtract all close solu\n")

    fi.write("\n")
    fi.write('print "beginning rerun %d"\n'%frame)
    fi.write("rerun ../../dump_dir/dump.lammpsdump-%d first %d last %d dump x y z\n" % (frame,frame,frame))
    fi.write('print "finishing rerun"\n')
    fi.close()

def Write_System(molecules, sysfile="system.in", subdir="rdf_test"):
    """This writes the includes into the simulation input file.

    This takes a base system input file (Iargs.sys) and writes them into a separate system input file
    for each molecule, so that lammps can easily calculate all of them.

    Args:
        nmols (int): Number of molecules
        sysfile (str): System input file with no run commands inside of it.
        subdir (str): Subdirectory where rdf calculation will take place

    """
    print("Writing System Files")
    for molid in molecules:
        with open("./%s/calc_dir/%s-%d"%(subdir,sysfile,molid),'w') as g:
            g.write("include ../include_header\n")
            g.write("log logs/log.%d.out\n"%molid)
            g.write("include ../include_dir/include.computes\n")
            g.write("include ../include_dir/include.groups-%d\n"%molid)
    
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
            f.write("   molid=$((fid*%d + i ))\n" % nmol_per_task)
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
    parser.add_argument('-step', default=1, type=int, help='Which step of the calc to run')
    parser.add_argument('-nmols',   default=343,                  type=int, help='Number of molecules in the system')
    parser.add_argument('-nconfigs',default=10000,              type=int, help='Number of configurations [default 10000]')
    parser.add_argument('-sys',     default="system.in",        type=str, help='System input file [default system.in]')
    parser.add_argument('-skip',    default=1,                  type=int, help='Frequency to use frames')
    parser.add_argument('-framesep',default=1000,               type=int, help='Dump Frequency (Excluding Skip)')
    parser.add_argument('-startframe',default=1000000,          type=int, help='Starting frame timestep number')
    parser.add_argument('-subdir',  default='rdf_deriv',        type=str, help='Subdirectory for calculation')
    Iargs = parser.parse_args()

    if Iargs.step == 1:
        Prep_Dir(Iargs.nmols, Iargs.subdir)
    elif Iargs.step == 2:
        molecules = np.arange(Iargs.nmols)+1
        LAMMPS_Fort_Main(molecules, nconfigs=Iargs.nconfigs, skip=Iargs.skip, framesep=Iargs.framesep,
                         startframe=Iargs.startframe, subdir=Iargs.subdir)
    elif Iargs.step == 3:
        Write_Task(nmols = Iargs.nmols, nmol_per_task = 1, hours = 24, queue_engine = "TORQUE", subdir = Iargs.subdir)
    else:
        raise ValueError("Error: Step %d is not yet supported" % Iargs.step)