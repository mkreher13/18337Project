#18.337 Project RandomRay Miriam Kreher

#Imports
from mpi4py import MPI
from Geometry import *
from time import time
startTime = time()

#Variables
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
name = MPI.Get_processor_name()

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
tmp = 1./4./np.pi

idx = range(rays)
myidx = idx[rank:rays:size]
print("Hello World! I am process %d of %d on %s.\n" % (rank, size, name))


#Initial functions
G = Geometry(pitch, rays)
#Initialize cross sections, source, etc
F = Flux(fuel_rad, pitch, nCells, nGrps, NEWk, rad, rays)

#Iterations 
while abs(k-NEWk)>1e-5 or m>1e-5:
# for i in range(0,1):

	#Updates for this iteration
	iteration = iteration+1
	# random.seed(1)
	k = NEWk
	F.deltapsi_storage = np.zeros([rays, nCells,nGrps])
	F.phibar = np.zeros([nCells,nGrps])
	G.vol_storage = np.zeros([rays,nCells])
	# G.vol = np.zeros(nCells)

	#Build rays
	for ray in myidx:
	# for ray in range(0,rays):
		#Count distance in moderator, fuel, total distance
		d_tot = 0
		#Trace rays
		G.ray_tracing(F, ray_length, deadzone, d_tot, nGrps, rad, ray)
	my_storage = np.sum(F.deltapsi_storage,axis=0)
	my_vol = np.sum(G.vol_storage,axis=0)
	# print("Hello World! I am process ",rank," My value is: ",np.sum(F.deltapsi_storage,axis=0))
	# print(np.sum(F.deltapsi_storage,axis=0))

	comm.Barrier()
	all_storage = comm.gather(my_storage, root=0)
	all_vol = comm.gather(my_vol, root=0)
	if rank == 0:
		deltapsi_storage = np.zeros([nCells,nGrps])
		volume = np.zeros(nCells)
		for storage in all_storage:
			deltapsi_storage = deltapsi_storage + storage 
		for vol in all_vol:
			volume = volume + vol
	else:
		deltapsi_storage = None
		volume = None
	deltapsi_storage = comm.bcast(deltapsi_storage, root=0)
	volume = comm.bcast(volume, root=0)


	#Add one-time terms to phibar
	for cell in range(0,nCells):
		F.phibar[cell,:] = F.phibar[cell,:]+4.*np.pi*deltapsi_storage[cell,:] #np.sum(F.deltapsi_storage,axis=0)[cell,:]
		F.phibar[cell,:] = F.phibar[cell,:]/volume[cell] #(G.vol[cell])
	F.phibar[:,:] = F.phibar[:,:]/F.tot[:,:]+F.q[:,:]*4*np.pi


	#Update eigenvalue k
	power = sum(F.oldphi[:,:]*F.nufiss[:,:])
	NEWpower = sum(F.phibar[:,:]*F.nufiss[:,:])
	NEWk = sum(NEWpower)/sum(power)*k
	
	#Normalize flux phibar
	F.phibar[:,:] = F.phibar[:,:]/np.linalg.norm(F.phibar)

	#Update source q
	nuf_rate = F.phibar[:,:]*F.nufiss[:,:]
	for cell in range(0,nCells):
		for g in range(0,nGrps):
			F.q[cell,g] = tmp/F.tot[cell,g]*(F.chi[cell,g]/NEWk*sum(nuf_rate[cell,:]
				)+sum(F.scat[:,g+nGrps*cell]*F.phibar[cell,:]))


	#Calculate convergence criteria
	m = np.ndarray.max(F.q[:,:]-F.oldq[:,:])
	F.oldq = copy.copy(F.q)
	F.oldphi = copy.copy(F.phibar)


#Output convergence criteria
if rank == 0:
	print('Source:') 
	print(F.q)
	print('k: %f' % (NEWk))

	metric = (time()-startTime) #/(F.seg_counter*nGrps)
	print('Execution metric [time] is %.8f seconds' %(metric))


# print('Number of iterations is %i' % (iteration))
# print(F.seg_counter)
# metric = (time()-startTime) #/(F.seg_counter*nGrps)
# print('Execution metric [time] is %.8f seconds' %(metric))

