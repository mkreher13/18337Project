#18.337 Project Miriam Kreher
#Class to create & operate on geometry and rays

from Flux import *
import numpy as np
import random

class Geometry():

	def __init__(self, Pitch, rays):

		random.seed(1)

		self.RandNums = np.zeros([rays,4])
		for i in range(0,rays):
			for j in range(0,4):
				self.RandNums[i,j] = random.random()
		
		#defining surfaces 
		HalfPitch = Pitch/2.
		self.xRight = HalfPitch
		self.xLeft = -HalfPitch
		self.yTop = HalfPitch
		self.yBottom = -HalfPitch
	

######################################################	

	def ray_tracing(self, F, max_length, deadzone, d_tot, nGrps, rad, ray):

		#Sample new starting points & angle
		StartX = self.RandNums[ray,0]*1.26-0.63 #random.random()*1.26-0.63
		StartY = self.RandNums[ray,1]*1.26-0.63 #random.random()*1.26-0.63
		theta = self.RandNums[ray,2]*2*np.pi #random.random()*2*np.pi
		u = self.RandNums[ray,3] #random.random()

		#Reset deltapsi and psi
		d_center_start = np.sqrt(StartX**2+StartY**2)
		if d_center_start <= rad[0]:
			cell=1
			F.newray(nGrps, cell)
		else:
			cell=0
			F.newray(nGrps, cell)

		#Build segments
		while d_tot < max_length:

			#D is list of distances from all boundaries
			#Distances added only if strictly positive & real
			D = []
			d5 = 0
			d6 = 0
			#d1 is distance from right boundary
			d1 = (self.xRight-StartX)/np.cos(theta)
			if d1 > 1e-14:
				D.append(d1)
			#d2 is distance from left boundary
			d2 = (self.xLeft-StartX)/np.cos(theta)
			if d2 > 1e-14:
				D.append(d2)
			#d3 is distance from top boundary
			d3 = (self.yTop-StartY)/np.sin(theta)
			if d3 > 1e-14:
				D.append(d3)
			#d4 is distance from bottom boundary
			d4 = (self.yBottom-StartY)/np.sin(theta)
			if d4 > 1e-14:
				D.append(d4)
			#d5 and d6 are distances from all circle intersections
			k = StartX*np.cos(theta)+StartY*np.sin(theta)
			#Iterate over spatially refined cirles within fuel
			for r in rad:
				c = StartX**2+StartY**2-r**2
				if k**2-c > 0:
					d5 = -k+np.sqrt(k**2-c)
					if d5 > 1e-14:
						D.append(d5)
					d6 = -k-np.sqrt(k**2-c)
					if d6 > 1e-14:
						D.append(d6)

			#Select shortest distance
			d = min(D)
			s = d/u
			d_tot = d_tot+s
			#Determine end points of segment
			EndX = StartX+d*np.cos(theta)
			EndY = StartY+d*np.sin(theta)

			#Determine which cell the segment is in
			#and calculate its flux contribution
			#and change angle if it hits cell boundary
			if d == d1 or d == d2:
				#Hits a right or left boundary
				theta = -theta+np.pi
				cell=0
				if d_tot > deadzone:
					F.contribute(cell,s,nGrps,ray)
					self.vol_storage[ray,cell] = self.vol_storage[ray,cell]+s
					# self.vol[cell] = self.vol[cell]+s
			elif d == d3 or d == d4:
				#Hits a top or bottom boundary
				theta = -theta
				cell=0
				if d_tot > deadzone:
					F.contribute(cell,s,nGrps,ray)
					self.vol_storage[ray,cell] = self.vol_storage[ray,cell]+s
					# self.vol[cell] = self.vol[cell]+s
			else:
				#Distance from the center for the start point of the ray and
				#distance from the center for the end point of the ray
				#to determine which cell the segment is in
				d_center_start = np.sqrt(StartX**2+StartY**2)
				d_center_end = np.sqrt(EndX**2+EndY**2)
				a = round(d_center_start,1)
				b = round(d_center_end,1)
				#Assume cell is moderator, and correct if 
				#location is determined to be in fuel
				cell = 0
				for i in range(0,len(rad)):
					if round(d_center_start,3) == rad[i]:
						if a == b:
							cell = i+1
						elif a < b:
							cell = i
						elif a > b:
							cell = i+1
				if d_tot > deadzone:
					F.contribute(cell,s,nGrps,ray)
					self.vol_storage[ray,cell] = self.vol_storage[ray,cell]+s
					# self.vol[cell] = self.vol[cell]+s


			StartX = EndX
			StartY = EndY
