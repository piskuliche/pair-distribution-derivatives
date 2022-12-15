#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt

"""
Copyright November 2022, Boston University

This is a class that stores information about RDFs for each molecule. It also contains class 
methods to weight the RDF by energies either on anindividual basis, or on a collective basis. 

For questions, contact Zeke Piskulich, piskulichz@gmail.com

"""

class RDFS:
    """Calculates and stores the RDF

    This class calculates and stores the RDF, and can calculate the derivative of the RDF 
    if it is provided with energies, calculated elsewhere.

    Attributes:
        dr (float): Bin width in A for RDF 
        rmax (float): Cutoff distance in A for the RDF 
        nbins (int): Number of RDF bins
        dim (int): Number of dimensions 
        rdfs (dict): Dictionary of lists, each dict element is a different molecule
        xrdf (array): Numpy array of x values, used for plotting

    """

    def __init__(self, dr=0.1, rmax=12, dim = 2):
        """Sets up the RDF class and builds attributes

        This method takes user input and begins the initialization of the RDFS class.
        This can be eventually extended, or split into smaller classes using class
        methods.

        Args:
            dr (float): Bin width in A for RDF [default=0.1]
            rmax (float): Cutoff distance in A for the RDF [default=12]
            dim (int): Number of dimensions [default=2]

        """
        self.dr     = dr
        self.rmax   = rmax
        self.nbins  = int(rmax/dr)
        self.dim    = dim
        self.rdfs   = {}
        self.xrdf   = np.arange(0,self.rmax,self.dr) + self.dr/2

    def calc_rdf(self, r, L, ids, init=False):
        """Function to calculate the RDF

        This is a class for calculating the rdf, given a set of positions. Note that self.dim chooses between a 2D or 3D RDF.

        Args:
            r (array_like): Coordinates of the atoms for which the RDF is being calculated.
            L (array_like): Box dimension array
            ids (array_like): Array of molecule ids to act as dictionary keys for self.rdfs

        """

        natoms = np.shape(r)[0]
        rlo    = np.arange(0,self.rmax,self.dr)
        rhi    = rlo + self.dr

        # Set h_id based on the dimensionality
        h_id = 1.0
        if self.dim == 2:   
            h_id = math.pi*natoms/(L[0]*L[1])*(rhi**2.-rlo**2.)
        elif self.dim == 3:
            h_id = 4/3*math.pi*natoms/(L[0]*L[1]*L[2])*(rhi**3. - rlo**3.)
        
        # Builds the dictionary if init is true
        if init == True:
            for i in range(natoms):
                assert ids[i] not in self.rdfs
                self.rdfs[ids[i]] = []

        # loops over atoms, twice. Right now, only same rdfs works. 
        for i in range(natoms):
            counts = np.zeros(self.nbins)
            for j in range(natoms):
                if i!=j:
                    dist=self._pbc_dist(r[i],r[j],L)
                    if dist < self.rmax:
                        counts[self._find_nearest(self.xrdf,dist)] += 1
            self.rdfs[ids[i]].append(counts/h_id)


    def report_rdfs(self, plot=False, average_over_configs=True):
        """Averages the RDF and plots it

        This function returns the radial distribution function averaged over all molecules and configurations in the 
        stored class. It also plots the RDF for your visual pleasure.

        Args:
            plot (bool): True, plots the RDF. False, does not plot. [default=False]
            average_over_configs (bool): True, returns the config averaged rdf, False, returns pre-config averaging.

        Returns:
            (tuple): RDF(s) and its Uncertainty
            
                * avrdf (array_like): average_over_configs=True: Averaged RDF as a function of distance. 
                    average_over_configs = False, array of rdfs as a function of distance shape (nconfigs,r)
                * err (array_like): 95% CI for the RDF based on nblocks

        """
        all_mols = []
        for key in self.rdfs:
            all_mols.append(self.rdfs[key])
        av_over_mol = np.average(all_mols,axis=0)
        av_total = np.average(av_over_mol,axis=0)
        err = self._Error(av_over_mol,nblocks=5)
        
        if plot == True:
            fig = plt.figure(dpi=300,figsize=(3,3))
            plt.plot(self.xrdf,av_total)
            plt.xlabel("r (Angstroms)")
            plt.ylabel("g(r)")
            plt.savefig("RDF_reported.png")

        if average_over_configs == True:
            return av_total, err
        else:
            return av_over_mol, err
    
    def total_ener_weight(self, totener, nblocks=1, average_energy = 0.0, plot=False, average_over_configs=True):
        """Weight RDF by total energy fluctuation.

        This implements the fluctuation theory approach to calculating a derivative and applies it to the
        RDF to obtain the negated derivative of the RDF. This does energy weighting by total energies, rather
        than for individual moelucles. If you want to do individual weighting, use mol_ener_weight instead.

        The feature, average_over_configs, sets whether the array_like return for the RDF derivative is
        averaged over configurations. If True, it is averaged over configurations. If False, it is not. This 
        can be useful for when you want to break the configuration up differently for block averaging over
        multiple class objects.

        Average energies can be overridden by the average_energy keyword. By default, the code calculates the average energy itself,
        but if desired, this can be supplied by the user instead through this keyword.

        Args:
            totener (array_like): Array of systemwide energies, must match number of configs.
            nblocks (int): Number of blocks for block averaging [default=1]
            plot (bool): True, plots the derivative. False, does not plot. [default=False]
            average_over_configs (bool): True, does config averaging. Flase, does not. [default=True]
            average_energy (float): Value of the average energy, if non-zero, will use this instead of the average energy
                calculated from totener. This is useful for cases like block averaging.

        Returns:
            (tuple): Contains the RDF output and the error.
                
                * avrdf (array_like): average_over_configs=True: Averaged dRDF as a function of distance. 
                    average_over_configs = False, array of drdfs as a function of distance shape (nconfigs,r)
                * err (array_like): 95% CI for the RDF based on nblocks

        """

        if average_energy != 0.0: average_energy = np.average(totener)

        all_mols=[]
        for key in self.rdfs:
            all_mols.append(self.rdfs[key])
        av_over_mol = np.average(all_mols,axis=0)
        de = totener - average_energy
        wrdf = np.multiply(de[:,None],av_over_mol)
        av_total = np.average(wrdf,axis=0)

        err = self._Error(av_over_mol,nblocks=nblocks)

        if plot == True:
            fig = plt.figure(dpi=300,figsize=(3,3))
            plt.plot(self.xrdf, wrdf)
            plt.xlabel("r (Angstroms)")
            plt.ylabel("g_H(r)")
            plt.savefig("Deriv_RDF_reported.png")
        
        if average_over_configs == True:
            return av_total, err
        else:
            return av_over_mol, err


    def mol_ener_weight(self, ener, nblocks=1, average_energy=0.0, plot=False, average_over_configs=True):
        """Weight RDF by molecular energy fluctuation. 

        This implements the fluctuation theory approach to calculating a derivative and applies it to
        the RDF to obtain the negated derivative of the RDF. This does energy weighting by molecule, rather
        than for the whole system. If you want to do the whole system, use total_ener_weight instead.

        The feature, average_over_configs, sets whether the array_like return for the RDF derivative is
        averaged over configurations. If True, it is averaged over configurations. If False, it is not. This 
        can be useful for when you want to break the configuration up differently for block averaging over
        multiple class objects.

        Average energies can be overridden by the average_energy keyword. By default, the code calculates the average energy itself,
        but if desired, this can be supplied by the user instead through this keyword.

        Args:
            ener (dict): Dictionary that stores energies for each molecule, for each configuration. 
                Dictionary keys must match self.rdfs, and must have same number of configurations.
            nblocks (int): Number of blocks for block averaging (default=1)
            plot (bool): True, plots the derivative. False, does not plot. [default=False]
            average_over_configs (bool): True, does config averaging. Flase, does not. [default=True]
            average_energy (float): Value of the average energy, if non-zero, will use this instead of the average energy
                calculated from ener. This is useful for cases like block averaging.
        
        Returns:
            (tuple): Contains the RDF derivative output and the error.
                
                * avrdf (array_like): average_over_configs=True: Averaged dRDF as a function of distance. 
                    average_over_configs = False, array of drdfs as a function of distance shape (nconfigs,r)
                * err (array_like): 95% CI for the RDF based on nblocks
                
        """
        if average_energy != 0.0: average_energy = np.average(list(ener.values()))

        all_mols = []
        for key in self.rdfs:
            de = ener[key] - average_energy
            wrdfs = np.multiply(de[:,None],self.rdfs[key])
            all_mols.append(wrdfs)

        av_over_mol = np.average(all_mols,axis=0)
        av_total = np.average(av_over_mol,axis=0)
        err = self._Error(av_over_mol,nblocks=nblocks)

        if plot == True:
            fig = plt.figure(dpi=300,figsize=(3,3))
            plt.plot(self.xrdf,av_total)
            plt.xlabel("r (Angstroms)")
            plt.ylabel("g_H(r)")
            plt.savefig("Deriv_RDF_reported.png")

        if average_over_configs == True:
            return av_total, err
        else:
            return av_over_mol, err
    
    def lightweight_save(self,subdivide=5000,subdir="rdf_classes"):
        """Splits the RDF class into smaller classes

        This class function takes the current class, and rebuilds it as a series of smaller classes. This is done
        by taking the total number of rdfs and dividing it into chunks of [subdivide] rdfs long. 
        
        For instance, if you have 50k rdfs, and subdivide is 5000, you get 10 class instances saved. 

        This lets you work with class objects, not in memory, which can be ideal. 

        Args:
            subdivide (int): Number of configurations per class instance [default=5000]
                Total number of configurations should be evenly divided by this number.
            subdir (str): Name of a directory to store the classes. This directory will be created
                if it does not already exist.

        """
        import pickle, os 
        # Make the subdirectory if it doesn't already exist.
        if not os.path.exists("./%s"%subdir):
            os.makedirs("./%s"%subdir)
        # Find the keys in the original class, and setup division
        sub_rdfs, keys = [], []
        for key in self.rdfs:
            keys.append(key)
        num_rdfs = np.shape(self.rdfs[keys[0]])[0]
        num_bins = int(num_rdfs/subdivide)
        # Build the New Classes
        for i in range(num_bins):
            start_idx, stop_idx=i*subdivide, (i+1)*subdivide
            sub_rdfs.append(RDFS(dr=self.dr,rmax=self.rmax,dim=self.dim))
            for key in keys:
                sub_rdfs[i].rdfs[key] = self.rdfs[key][start_idx:stop_idx]
        # Save the new class files        
        count = 0
        for sub in sub_rdfs:
            pickle.dump(sub,open("%s/rdf_data_%d.pckl"%(subdir,count),'wb'))
            count += 1
        print("Lightweight Save Complete")


    def _pbc_dist(self, r1, r2, L):
        """Calculates minimum image distance

        This function takes two vectors and calculates the minimum image distance by taking
        periodic boundary conditions into account
        
        Args:
            r1, r2 (array_like): two sets of coordinates of shape (3,)
            L (array_like): Box dimension array

        Returns: 
            float: Minimum image distance

        """
        dr, drsq = np.zeros(self.dim), 0.0
        for i in range(self.dim):
            dr[i] =  r1[i] - r2[i]
            dr[i] -= L[i]*round(dr[i]/L[i])
            drsq  += dr[i]**2.
        dist = np.sqrt(drsq)
        return dist

    def _find_nearest(self, array, value):
        """Find nearest location to the current value

        Args:
            array (array_like): Array to find value in
            value (float): Value to be located

        Returns: 
            int: Index in array for which value is located
            
        """
        array   = np.asarray(array)
        idx     = (np.abs(array-value)).argmin()
        return idx

    def _Error(self, data, nblocks=5):
        """Calculates the block average of a 2D array

        This function takes a 2D array (configs, M), and splits configs into
        nblocks blocks. 

        Args:
            data (array_like): 2D Array of data values of shape(configs,M)
            nblocks (int): Number of blocks for block averaging

        Returns:
            array_like: 95% CI of shape(M)

        """
        from scipy import stats
        if nblocks > 1:
            t_val=stats.t.ppf(0.975,nblocks-1)/np.sqrt(nblocks)
            blocks = np.average(np.array_split(data,nblocks),axis=1)
            err = np.std(blocks,axis=0)*t_val
        else:
            err = np.zeros(np.shape(data)[-1])
        return err
        
