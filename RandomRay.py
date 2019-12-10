#18.337 Project RandomRay Miriam Kreher

#Imports
from Geometry import *
from time import time
startTime = time()

#Variables
rays = 200 #200
ray_length = 75 #[cm]
deadzone = 25 #[cm]
pitch = 1.26
fuel_rad = 0.4
rad = np.array([0.4, 0.3, 0.2, 0.1])
nCells = len(rad)+1
nGrps = 2
NEWk = 1.0
k = 2.0
m = 2.0
iteration = 0

#Initial functions
G = Geometry(pitch)
#Initialize cross sections, source, etc
F = Flux(fuel_rad, pitch, nCells, nGrps, NEWk, rad)

#Iterations 
while abs(k-NEWk)>1e-5 or m>1e-5:

	#Updates for this iteration
	iteration = iteration+1
	random.seed(1)
	k = NEWk
	F.phibar = np.zeros([nCells,nGrps])
	G.vol = np.zeros(nCells)
	G.s_tot = 0

	#Build rays
	for ray in range(0,rays):
		#Count distance in moderator, fuel, total distance
		d_tot = 0
		#Trace rays
		G.ray_tracing(F, ray_length, deadzone, d_tot, nGrps, rad)

	#Add one-time terms to phibar
	for cell in range(0,nCells):
		F.phibar[cell,:] = F.phibar[cell,:]/(G.vol[cell])
	F.phibar[:,:] = F.phibar[:,:]/F.tot[:,:]+F.q[:,:]*4*np.pi

	#Update eigenvalue k
	power = sum(F.oldphi[:,:]*F.nufiss[:,:])
	NEWpower = sum(F.phibar[:,:]*F.nufiss[:,:])
	NEWk = sum(NEWpower)/sum(power)*k
	
	#Normalize flux phibar
	F.phibar[:,:] = F.phibar[:,:]/np.linalg.norm(F.phibar)

	#Update source q
	nuf_rate = F.phibar[:,:]*F.nufiss[:,:]
	tmp = 1./4./np.pi
	# scatsrc = np.zeros([nCells,nGrps])
	for cell in range(0,nCells):
		for g in range(0,nGrps):
			# scatsrc[cell,g] = sum(F.scat[:,g+nGrps*cell]*F.phibar[cell,:])
			F.q[cell,g] = tmp/F.tot[cell,g]*(F.chi[cell,g]/NEWk*sum(nuf_rate[cell,:]
				)+sum(F.scat[:,g+nGrps*cell]*F.phibar[cell,:]))
	# for cell in range(0,nCells):
	# 	for g in range(0,nGrps):
	# 		F.q[cell,g] = tmp/F.tot[cell,g]*(F.chi[cell,g]/NEWk*sum(nuf_rate[cell,:]
	# 			)+scatsrc[cell,g])

	#Calculate convergence criteria
	m = np.ndarray.max(F.q[:,:]-F.oldq[:,:])
	F.oldq = copy.copy(F.q)
	F.oldphi = copy.copy(F.phibar)

	#Output convergence criteria
	# print('Source:') 
	# print(F.q)
	# print('k: %f' % (NEWk))

# print('Number of iterations is %i' % (iteration))
metric = (time()-startTime)/(F.seg_counter*nGrps)
print('Execution metric [time/seg*nGrps] is %.8f microseconds' %(metric*1000000))

