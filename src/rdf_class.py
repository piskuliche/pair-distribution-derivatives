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
        self.dr     = dr
        self.rmax   = rmax
        self.nbins  = int(rmax/dr)
        self.dim    = dim
        self.rdfs   = {}

    def _calc_rdf(self, r, ts,ids,init=False):
        natoms = np.shape(r)[0]
        L      = ts.dimensions[:self.dim]

        h_real = np.zeros(self.nbins)
        self.xrdf   = np.arange(0,self.rmax,self.dr) + self.dr/2
        rlo    = np.arange(0,self.rmax,self.dr)
        rhi    = rlo + self.dr

        h_id = math.pi*natoms/(L[0]*L[1])*(rhi**2.-rlo**2.)
        
        if init == True:
            for i in range(natoms):
                assert ids[i] not in self.rdfs
                self.rdfs[ids[i]] = []
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
        all_mols = []
        for key in self.rdfs:
            all_mols.append(np.average(self.rdfs[key],axis=0))
        av=np.average(all_mols,axis=0)
        fig = plt.figure()
        plt.plot(self.xrdf,av)
        plt.show()
        return av

    def ener_weight(self,ener):
        all_mols = []
        for key in self.rdfs:
            eav = np.average(ener[key])
            de = ener[key] - eav
            wrdfs = np.multiply(de[:,None],self.rdfs[key])
            all_mols.append(np.average(wrdfs,axis=0))
        av = np.average(all_mols,axis=0)
        return av
        

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
