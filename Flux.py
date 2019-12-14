#18.337 Project Miriam Kreher
#Class to calculate flux contributions

import numpy as np
import copy

class Flux():

	def __init__(self, fuel_rad, P, nCells, nGrps, NEWk, rad, rays):

		#Initialize and allow for arbitrary amount of spatial refinement
		# self.seg_counter = 0
		# self.deltapsi_storage = {}



		# for ray in range(0,rays):
		# 	for cell in range(0,nCells):
		# 		for g in range(0,nGrps):
		# 			self.d = {
		# 			    ray: {
		# 			        cell: {
		# 			            g: {}
		# 			        }
		# 			    }
		# 			}
		# 			self.deltapsi_storage.update(self.d)

		# # print(self.deltapsi_storage)
		# self.deltapsi_storage[0][1][0] = 1.0
		# # self.deltapsi_storage[0,0,1] = 2.0
		# print(self.deltapsi_storage)

		# self.deltapsi_storage = np.zeros([rays, nCells,nGrps])

		self.vol = np.zeros(nCells)
		self.vol[0] = P**2-np.pi*rad[0]**2
		self.vol[1::] = np.pi*rad[:]*rad[:]

		self.abs = np.zeros([nCells,nGrps])
		self.abs[:,:] = [0.01,0.2]
		self.abs[0,:] = [0.01,0.02]

		self.tot = np.zeros([nCells,nGrps])
		self.tot[:,:] = [0.37588,0.634259]
		self.tot[0,:] = [0.6255565,1.8116]

		self.phibar = np.zeros([nCells,nGrps])
		self.phibar[:] = 1.
		self.oldphi = copy.copy(self.phibar)

		self.nufiss = np.zeros([nCells,nGrps])
		self.nufiss[:,:] = [0.013457,0.31337]
		self.nufiss[0,:] = [0.0,0.0]

		self.scat = np.zeros([nGrps,nGrps*nCells])
		self.scat[0,0] = 0.563958
		self.scat[0,1] = 0.062163
		self.scat[0,2::2] = 0.3643677
		self.scat[0,3::2] = 0.0012458
		self.scat[1,1] = 1.7964575
		self.scat[1,3::2] = 0.42842

		self.q = np.zeros([nCells,nGrps])

		self.chi = np.zeros([nCells,nGrps])
		self.chi[:,:] = [1.0,0.0]
		self.chi[0,:] = [0.0,0.0]


		nuf_rate = self.phibar[:,:]*self.nufiss[:,:]
		tmp = 1./4./np.pi
		for cell in range(0,nCells):
			for g in range(0,nGrps):
				self.q[cell,g] = tmp/self.tot[cell,g]*(self.chi[cell,g]/NEWk*sum(nuf_rate[cell,:]
					)+sum(self.scat[:,g+nGrps*cell]*self.phibar[cell,:]))
		self.oldq = copy.copy(self.q)
		

######################################################

	def newray(self, nGrps, cell):
		self.deltapsi = np.zeros(nGrps)
		self.psi = copy.copy(self.q[cell])

######################################################

	def contribute(self, cell, s, nGrps, ray):

		# self.seg_counter = self.seg_counter+1

		for g in range(0, nGrps):
			self.deltapsi[g] = (self.psi[g]-self.q[cell,g])*(
				1-np.exp(-self.tot[cell,g]*s))

			self.deltapsi_storage[ray,cell,g] = (
				self.deltapsi_storage[ray][cell][g] + self.deltapsi[g])

			# self.phibar[cell,g] = self.phibar[cell,g]+4*np.pi*self.deltapsi[g]
			self.psi[g] = self.psi[g]-self.deltapsi[g]
