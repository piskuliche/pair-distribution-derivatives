#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt

class Energies:
    """Class to store energies for energy weighting

    Attributes:
        mol_ener (dictionary): Dictionary of energy values
        tot_ener (dictionary): Dictionary of system-wide energies 
        nmols (int): Number of molecules

    """
    def __init__(self,nmols):
        """Builds the basics for the class

        Args:
            nmols (int): Number of molecules

        """
        self.mol_ener = {}
        self.tot_ener = {}
        self.nmols = nmols

    def add_mol_ener(self, compstr='ss', col=1, fname='pair.compute', program='LAMMPS'):
        """Add energies from a file to the self.mol_ener dictionary

        This function reads in energies from the file and gets it into the molecular energy container such that it
        could be used for later weighting. 

        Args: 
            compstr (str): String name for the energy component
            col (int): Integer number for column from which energy is read.
            fname (str): Base name of energy file.
            program (str): Which MD program generated the files.

        """
        if compstr in self.mol_ener:
            raise ValueError("%s already exists within class."%compstr)
        self.mol_ener[compstr]={}
        for i in range(self.nmols):
            self.mol_ener[compstr][i+1]=np.genfromtxt("ener_dir/%s%d"%(fname,i),unpack=True,usecols=col,skip_header=2)
    
    def pull_lammps(self, fname='log.production'):
        """Pulls energy data from LAMMPS

        Pulls energy for all columns in the LAMMPS log file and writes them
        as a file for each one. This intelligently parses the lammps logfile,
        but must be a full logfile with only a single run command, otherwise breaks

        Args:
            fname (str): Name of the logfile.

        """
        data={}
        with open(fname, 'r') as f:
            lines=f.readlines()
            flag=0
            keys=[]
            for line in lines:
                if "Loop" in line:
                    flag=0
                    print("stop")
                if flag == 1:
                    for key in keys:
                        data[key].append(float(line.strip().split()[loc[key]]))
                if "Step Time" in line:
                    flag=1
                    data={}
                    loc={}
                    keys=line.strip().split()
                    count = 0
                    for key in keys:
                        data[key]=[]
                        loc[key]=count
                        count+=1
                    print("start")

        for key in data:
            np.savetxt("%s_init.out" % key, np.c_[data[key]])
            self.tot_ener[key] = data[key]
            if key == "Volume":
                L = np.array(data[key])**(1./3.)
                np.savetxt("L.dat", np.c_[L])

    def gen_default(self, fname='pair.compute',program='LAMMPS'):
        """Generates a default set of energies

        This function assumes that you want to split your energies by solute, close and far, and the intermingling of 
        these groups, rather than something more complex. It reads all the files, and generates the dictionary with them in it.

        Args:
            fname (str): Base name of energy file.
            program (str): Which MD program generated the files.

        """
        default_types = ['ss','cc','ff','sc','sf','cf']
        for i, key in enumerate(default_types):
            print("Reading energies for %s"%key)
            self.add_mol_ener(compstr=default_types[i], col=i+1, fname=fname, program=program)

    def plot_ener_dist(self,compstr='ss', _nbins=100, _dpi=300, _figsize=(3,3), savefile=False):
        """Plots the energy distributions

        This function plots the energy distributions using matplotlib.

        Args:
            compstr (str): String name for the energy component [default='ss']
            _nbins (int): Number of histogram bins [default=100]
            _dpi (int): DPI for plot [default=300]
            _figsize (tuple): Figure size dimensions [default=(3,3)]
            savefile (bool): Save to a file rather than print to screen [default=False]


        """
        e_list = []
        for key in self.mol_ener[compstr]:
            e_list.append(self.mol_ener[compstr][key])
        mol_av = np.average(e_list,axis=0)
        eav = np.average(e_list)

        de = mol_av - eav

        fig = plt.figure(dpi=_dpi, figsize=_figsize)
        plt.hist(de, bins=_nbins,density=True)
        plt.xlabel("$dE$ (kcal/mol)")
        plt.ylabel("P(dE)")
        plt.tight_layout
        if savefile == False:
            plt.show()
        else:
            plt.savefig("%s_distribution.png"%compstr)

    def values(self, compstr='ss', _start=0, _stop=-1, _skip=1):
        """Outputs the energy values, and splices them any which way

        This is just a simple way of pulling different energies (and subdividing them) into more
        manageable chunks. This is useful when the RDF class has been subdivided into smaller 
        class objects, while retaining a majority of the energies in one location.

        Args:
            compstr (str): String name for the energy component [default 'ss']
            _start (int): index to start splice [default=0]
            _stop (int): index to stop splice [default=-1]
            _skip (int): Skip between indices [default=1]

        Returns:
            dictionary: Energy values in a dictionary of arrays.

        """
        e_vals = {}
        for key in self.mol_ener[compstr]:
            e_vals[key]=self.mol_ener[compstr][key][_start:_stop:_skip]
        return e_vals

    def stats(self,compstr='ss', _start=0, _stop=-1, _skip=1):
        """Calculates and reports averages and statistics

        This operates statistics on a particular energy type.

        Args:
            compstr (str): String name for the energy component [default 'ss']
            _start (int): index to start splice [default=0]
            _stop (int): index to stop splice [default=-1]
            _skip (int): Skip between indices [default=1]
            
        """
        e_list = []
        for key in self.mol_ener[compstr]:
            e_list.append(self.mol_ener[compstr][key][_start:_stop:_skip])
        eav = np.average(e_list)
        estd = np.std(e_list)
        molav = np.average(e_list,axis=0)
        de = molav-eav
        deav = np.average(de)
        destd = np.std(de)

        print("--STATS--")
        print("Energy Type: %s" % compstr)
        print("<E>  %10.5f %10.5f" % (eav,estd))
        print("<dE> %10.5f %10.5f" % (deav,destd))
        print("--Done--")











