import numpy as np 
import matplotlib.pyplot as plt
import os.path 
from tqdm import tqdm
import time
from numba import jit
save_to_file=False
plot_histograms=False


def initialValues(file='parameters.txt'):
	global n,m,e,R,F,L,a,T0,tau,So,Sd,Sout,Sxyz,N,k
	initialVal=[]
	current_dir=os.path.dirname(__file__)
	file_path=os.path.join(current_dir,file)
	print(file_path)
	with open(file_path, 'r') as rFile: # the best way to open file which handles the exceptions
		for line in rFile: # for each line in file
			for val in line.split(): # split line
				try:
					initialVal.append(float(val))#float() returns exception when no valid number could be parsed
				except ValueError:
					pass
		initialVal[0]=int(initialVal[0])#number of particles is integer	
		n,m,e,R,F,L,a,T0,tau,So,Sd,Sout,Sxyz=initialVal 
		So=int(So)
		Sd=int(Sd)
		Sout=int(Sout)
		Sxyz=int(Sxyz)
		#L=1.22*a*(n-1)
		#n- number of atoms in one direction, m-atom mass, atom charge, R,f - elesticity coeficient for calculating force,L- diamater of force sphere,a- lattice constant, T0- temperature, tau,...
		N=n**3 # absolute number of atoms
		k=8.31e-3 # boltzman constant
		
def createInitialLattice(save=save_to_file):
	for ii0 in range(0,n):
		for ii1 in range(0,n):
			for ii2 in range(0,n):
				ii=ii0+ii1*n+ii2*n**2
				rtemp=(ii0-(n-1)/2)*B0+(ii1-(n-1)/2)*B1+(ii2-(n-1)/2)*B2
				#print(rtemp)
				#print('xd')
				r[ii,:]=rtemp# chooses i-th raw and changes all of the columns
				#print(r[ii])
	# write file in XYZ format to check in JMol if everything is ok
	if(save):
		data_file=open("r.xyz","w")
		data_file.write(f'{N}'+'\n')
		data_file.write('\n')
		for x,y,z in r:	
			data_file.write('Ar'+'\t'+f'{x}\t{y}\t{z}'+'\n')
		
def setInitialMomentum():
	#losowanie E kin
	Ekin=np.array([[0]*3 for ii in range (N)]).astype(np.float)
	for jj in range(0,N):
		for ww in range(0,3):
			rand=np.random.random_sample()
			if(rand==0):
				rand==1
			Ekin[jj][ww]=-1/2*k*T0*np.log(rand)
		plusMinus=np.array([np.random.choice([-1,1])]*3)
		#print(plusMinus)
		p[jj,:]=np.array([plusMinus[0]*np.sqrt(2*m*Ekin[jj][0]),plusMinus[1]*np.sqrt(2*m*Ekin[jj][1]),plusMinus[2]*np.sqrt(2*m*Ekin[jj][2])])
		#p[jj,:]=np.array([np.sqrt(2*m*Ekin[jj][nn]) if plusMinus[nn]>=0.5 else -1*np.sqrt(2*m*Ekin[jj][nn] for nn in range(0,3))])
	PP=np.array([0]*3).astype(np.float)
	PP=sum(p)
	PP=np.array([0]*3).astype(np.float)
	for jj in range(0,N):
		PP+=p[jj]
	print(PP)
	for jj in range(0,N):
		p[jj]=p[jj]-1/N*PP
	

def calculateFp_axis(r_temp,ax):
	Fp=np.array([[0]*N for ii in range (N)]).astype(np.float)
	for ii,ri in enumerate(r_temp):
		for jj,rj in enumerate(r_temp):
			if(jj>ii):#tutaj zmienilem warunek
				rij=np.linalg.norm(ri-rj)
				Fp[ii][jj]=12*e*((R/rij)**12-(R/rij)**6)*(ri[ax]-rj[ax])/rij**2
				Fp[jj][ii]=-Fp[ii][jj]
	return Fp

@jit(nopython=True, nogil=True)	
def calculateFp(r_temp):
	Fp=np.zeros((N,N,3))
	for ii, ri in enumerate(r_temp):
		for jj,rj in enumerate(r_temp[:ii]):# odrazu ii!=jj 
		#zmienilem
			#if(jj>ii):#tutaj zmienilem warunek
			if(jj!=ii):		
				rij=np.linalg.norm(ri-rj)
				Fp[ii][jj]=12*e*(np.power((R/rij),12)-np.power((R/rij),6))*(ri-rj)/np.power(rij,2)
				Fp[jj][ii]=-Fp[ii][jj]
	return np.sum(Fp,axis=1)

	
@jit(nopython=True, nogil=True)
def calculateFsLength(r_temp):
	Fsnorm=np.zeros(N)
	for ii,rr in enumerate(r_temp):
		ri=np.linalg.norm(rr)
		if( ri>L):
			Fsnorm[ii]=np.linalg.norm(1/ri*F*(L-ri)*rr)		
	return Fsnorm
	
	
def calculateVp(r_temp):
	vp=np.array([[0]*N for ii in range (N)]).astype(np.float)
	for ii,ri in enumerate(r_temp):
		for jj,rj in enumerate(r_temp):
			if(jj>ii):
				rij=np.linalg.norm(ri-rj)
				vp[ii][jj]=e*((R/rij)**12-2*(R/rij)**6)
	return vp
		
def calculateVs(r_temp):
	vs=np.array([0]*N).astype(np.float)
	for ii, ri in enumerate(r_temp):
		if(np.linalg.norm(ri)>L):
			vs[ii]=(1/2*F)*((L-np.linalg.norm(ri))**2)
			
	return vs
@jit(nopython=True, nogil=True)
def calculateForceAndPotential(r_temp):
	Vp=0
	Vs=0
	Fp=np.zeros((N,N,3))
	Fs=np.zeros((N,3))
	for ii, ri in enumerate(r_temp):
		r_abs=np.linalg.norm(ri)
		if( r_abs>L):
			Fs[ii,:]=1/r_abs*F*(L-r_abs)*ri
			Vs+=(1/2*F)*((L-np.linalg.norm(ri))**2)
		for jj,rj in enumerate(r_temp[:ii]):
			if(jj!=ii):		
				rij=np.linalg.norm(ri-rj)
				Fp[ii][jj]=12*e*(np.power((R/rij),12)-np.power((R/rij),6))*(ri-rj)/np.power(rij,2)
				Vp+=e*(np.power((R/rij),12)-2*np.power((R/rij),6))
				Fp[jj][ii]=-Fp[ii][jj]
	return Fs+np.sum(Fp,axis=1),Vs+Vp

def calculateWallPressure(r_temp):

	return 1/(4*np.pi*L**2)*calculateFsLength(r_temp).sum()

def calculateHamiltonian(p_temp,v_temp):

	return  calculateEk(p_temp)+v_temp
	
def calculateTemperature(p_temp):

	return 2/(3*N*k)*calculateEk(p_temp)
	
def calculateEk(p_temp):
	return ((np.power((np.linalg.norm(p_temp,axis=1)),2))).sum()/(2*m)

	
start_time=time.time()
initialValues() # loading initial constant values from file
#creating np array of dimensions 3xN
r=np.array([[0]*3 for ii in range (N)]).astype(np.float)
#print(r)
p=np.array([[0]*3 for ii in range (N)]).astype(np.float)
#positions of atoms in crystal corners
B0 = np.array([a, 0, 0])
B1 = np.array([a/2, a*np.sqrt(3)/2, 0])
B2 = np.array([a/2, a*np.sqrt(3)/6, a*np.sqrt(2/3)])
np.set_printoptions(precision=3)
#print(B0,B1,B2,sep=' ')
createInitialLattice() # preparing initial crystal lattice
setInitialMomentum() #setting initial momentum 
if(plot_histograms):
	# ploting histograms
	px=p[:,0]
	py=p[:,1]
	pz=p[:,2]
	fig, axs = plt.subplots(1,3)
	axs[0].set_title('px')
	axs[1].set_title('py')
	axs[2].set_title('pz')
	axs[0].hist(px,bins=20,label='px',color='r')
	axs[1].hist(py,bins=20,label='py',color='g')
	axs[2].hist(pz,bins=20,label='pz',color='b')
	plt.show()
	
	


force,v=calculateForceAndPotential(r)
print('Forces:')
print(force)
#calculate pressure exerted on sphere walls
print("V:\t"+f'{v}')
pressure=calculateWallPressure(r)
print("Pressure:\t"+f'{pressure}')


hamil=calculateHamiltonian(p,v)
print("Hamiltonian:\t"+f'{hamil}')
temperature=calculateTemperature(p)
print("Temperature:\t"+f'{temperature}')
#calculate momentum and position as function of t
#calculate momentum and position as function of t
#p_t=np.array([[0]*3 for ii in range (N)]).astype(np.float)
#p_t.append(p)
p_t=p
#r_t=np.array([[0]*3 for ii in range (N)]).astype(np.float)
r_t=r
#print("r(t)")
#print(r_t)
#print('p(t)')
#print(p_t)
f_tmp=force
file1=open("tHVEkTP.txt",'w')
file2=open("avs.txt",'w')
file3=open("avgTPH.txt",'w')
T_sum=0
P_sum=0
H_sum=0
t=0
for_start_time=time.time()
for jj in tqdm(range(0,So+Sd)):
	t=(jj+1)*tau
	p_tmp=p_t+1/2*f_tmp*tau
	r_t=(r_t+1/m*p_tmp*tau)
	#v_tmp=calculateVp(r_t)+calculateVs(r_t)
	f_tmp,v_tmp=calculateForceAndPotential(r_t)
	P_tmp=calculateWallPressure(r_t)
	#P_tmp=0
	p_t=(p_tmp+1/2*f_tmp*tau)
	T_tmp=calculateTemperature(p_t)
	H_tmp=calculateHamiltonian(p_t,v_tmp)
	
	
	if jj%Sout==0:
		#s1_start_time=time.time()
		ek_tmp=calculateEk(p_t)
		file1.write(f"{t}\t{H_tmp}\t{v_tmp}\t{ek_tmp}\t{T_tmp}\t{P_tmp}\n")
		#s1_end_time=time.time()
		#print('ts1'+f'{s1_start_time-s1_end_time}')
	if jj%Sxyz==0:
		#s2_start_time=time.time()
		file2.write(f'{N}'+'\n')
		file2.write('\n')
		for ii in range(0,N):
			x_tmp,y_tmp,z_tmp=r_t[ii][0],r_t[ii][1],r_t[ii][2]
			
			#file2.write(f"Ar\t{x_tmp:.2E}\t{y_tmp:.2E}\t{z_tmp:.2E}\t{ek_tmp:.2E}\n")
			file2.write(f"Ar\t{x_tmp}\t{y_tmp}\t{z_tmp}\n")
		#s2_end_time=time.time()
		#print('ts2'+f'{s2_start_time-s2_end_time}')
	
	if jj>=So:
		T_sum+=T_tmp
		P_sum+=P_tmp
		H_sum+=H_tmp
	
print('for loop took: ',for_start_time-time.time())

avgT=T_sum/Sd
avgP=P_sum/Sd
avgH=H_sum/Sd
sim_time=time.time()-start_time
file3.write(f"T_mean=\t{avgT:.2E}\tP_mean=\t{avgP:.2E}\tH_mean=\t{avgH:.2E}\n Simulation time: {sim_time}")		
file1.close()
file2.close()
file3.close()	
#print("r(t)")
#print(r_t)
#print('p(t)')
#print(p_t)