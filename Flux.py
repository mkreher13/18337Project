#18.337 Project Miriam Kreher
#Class to calculate flux contributions

import numpy as np
import copy

class Flux():

	def __init__(self, fuel_rad, P, nCells, nGrps, NEWk, rad):

		#Initialize and allow for arbitrary amount of spatial refinement
		self.phibar = []
		self.q = []
		self.tot = []
		self.abs = []
		self.chi = []
		self.nufiss = []
		self.vol = []
		self.scat = []
		self.phibar.append([1.,1.])
		self.q.append([0.,0.])
		self.tot.append([0.6255565,1.8116])
		self.abs.append([0.01,0.02])
		self.chi.append([0.0,0.0])
		self.nufiss.append([0.0,0.0])
		self.scat.append([0.563958,0.062163])
		self.scat.append([0.,1.7964575])
		for cell in range(0,nCells-1):
			self.phibar.append([1.,1.])
			self.q.append([0.,0.])
			self.tot.append([0.37588,0.634259])
			self.abs.append([0.01,0.2])
			self.chi.append([1.0,0.0])
			self.nufiss.append([0.013457,0.31337])
			for grp in range(0,nGrps):
				if grp == 0:
					self.scat[grp].append(0.3643677)
					self.scat[grp].append(0.0012458)
				else:
					self.scat[grp].append(0.)
					self.scat[grp].append(0.42842)
		for r in rad:
			if r == fuel_rad:
				self.vol.append(P**2-np.pi*r**2)
				self.vol.append(np.pi*r**2)
			else:
				self.vol.append(np.pi*r**2)
		self.phibar = np.array(self.phibar)
		self.q = np.array(self.q)
		self.tot = np.array(self.tot)
		self.abs = np.array(self.abs)
		self.chi = np.array(self.chi)
		self.nufiss = np.array(self.nufiss)
		self.vol = np.array(self.vol)
		self.scat = np.array(self.scat)

		nuf_rate = self.phibar[:,:]*self.nufiss[:,:]
		scatsrc = np.zeros([nCells,nGrps])
		for cell in range(0,nCells):
			for g in range(0,nGrps):
				scatsrc[cell,g] = sum(self.scat[:,g+nGrps*cell]*self.phibar[cell,:])
		for cell in range(0,nCells):
			for g in range(0,nGrps):
				self.q[cell,g] = 1./4./np.pi/self.tot[cell,g]*(self.chi[cell,g]/NEWk*sum(nuf_rate[cell,:]
					)+scatsrc[cell,g])
		self.oldq = copy.copy(self.q)
		self.oldphi = copy.copy(self.phibar)

		self.seg_counter = 0

######################################################

	def newray(self, nGrps, cell):
		self.deltapsi = np.zeros(nGrps)
		self.psi = copy.copy(self.q[cell])				

######################################################

	def contribute(self, cell, s, nGrps):

		self.seg_counter = self.seg_counter+1

		for g in range(0, nGrps):
			self.deltapsi[g] = (self.psi[g]-self.q[cell,g])*(
				1-np.exp(-self.tot[cell,g]*s))
			self.phibar[cell,g] = self.phibar[cell,g]+4*np.pi*self.deltapsi[g]
			self.psi[g] = self.psi[g]-self.deltapsi[g]
