#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt

"""
Copyright November 2022, Boston University

This is a class that stores information about RDFs for each molecule. It also contains class methods to weight the RDF by energies either on anindividual basis, or on a collective basis. 

For questions, contact Zeke Piskulich, piskulichz@gmail.com
"""

class RDFS:
    def __init__(self, dr=0.1, rmax=12, dim = 2):
        """
        Initializes data structures for the rdf
        """
        self.dr     = dr
        self.rmax   = rmax
        self.nbins  = int(rmax/dr)
        self.dim    = dim
        self.rdfs   = {}
        self.xrdf   = np.arange(0,self.rmax,self.dr) + self.dr/2

    def _calc_rdf(self, r, ts,ids,init=False):
        """
        This is a class for calculating the rdf. To do this, you provide it with positions (r),
        the timestep (ts), and a list of molecule ids (ids). 
        """

        natoms = np.shape(r)[0]
        L      = ts.dimensions[:self.dim]

        h_real = np.zeros(self.nbins)
        rlo    = np.arange(0,self.rmax,self.dr)
        rhi    = rlo + self.dr

        h_id = math.pi*natoms/(L[0]*L[1])*(rhi**2.-rlo**2.)
        
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
        return

    def report_rdfs(self):
        """
        This function returns the radial distribution function averaged over all molecules and configurations in the 
        stored class. It also plots the RDF for your visual pleasure.
        """
        all_mols = []
        for key in self.rdfs:
            all_mols.append(np.average(self.rdfs[key],axis=0))
        av=np.average(all_mols,axis=0)
        fig = plt.figure()
        plt.plot(self.xrdf,av)
        plt.show()
        return av

    def ener_weight(self,ener):
        """
        This function iterates through molecules and weights it by its energy fluctuation. 

        This outputs a single array the same length as the rdf, which corresponds to the weighted RDF 
        (and thus, the negated derivative)
        """
        all_mols = []
        for key in self.rdfs:
            eav = np.average(ener[key])
            de = ener[key] - eav
            wrdfs = np.multiply(de[:,None],self.rdfs[key])
            all_mols.append(np.average(wrdfs,axis=0))
        av = np.average(all_mols,axis=0)
        return av
    
    def lightweight_save(self,subdivide=5000,subdir="rdf_classes"):
        """
        This class function takes the current class, and rebuilds it as a series of smaller classes. This is done
        by taking the total number of rdfs and dividing it into chunks of [subdivide] rdfs long. 
        
        For instance, if you have 50k rdfs, and subdivide is 5000, you get 10 class instances saved. 

        This lets you work with class objects, not in memory, which can be ideal. 
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
        return

    def _pbc_dist(self, r1, r2, L):
        """
        This function takes two vectors and calculates the minimum image distance
            Returns: minimum image distance
        """
        dr, drsq = np.zeros(self.dim), 0.0
        for i in range(self.dim):
            dr[i] =  r1[i] - r2[i]
            dr[i] -= L[i]*round(dr[i]/L[i])
            drsq  += dr[i]**2.
        dist = np.sqrt(drsq)
        return dist

    def _find_nearest(self,array, value):
        """
        Finds the nearest location in array to value.
            Returns: index of location
        """
        array   = np.asarray(array)
        idx     = (np.abs(array-value)).argmin()
        return idx
