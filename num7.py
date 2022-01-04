import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm
def func(x):
	return np.sin(x*np.pi)


def EulerExplicit(dt,tmax,dx,xmax,D):
	r=D*dt/(dx*dx)
	n=int(xmax/dx)
	m=int(tmax/dt)
	U=np.zeros([n,m])
	for ii in range(0,m,1):
		U[0][ii]=0
		U[n-1][ii]=0
	t=np.arange(0,tmax,dt)
	x=np.arange(0,xmax,dx)
	for ii in range(1,n-1):
		x[ii]=ii*dx;
		U[ii][0]=func(x[ii])
	for jj in range(0,m-1):
		t[jj]=jj*dt
		for ii in range(1,n-1):
			U[ii][jj+1]=(r*U[ii-1][jj]+(1-2*r)*U[ii][jj]+r*U[ii+1][jj])	
	return x,t,U	
		
	
tmax=1
xmax=1
dt=tmax/2500
dx=xmax/30
r=dt/(dx*dx)
if(r<=1/2):
	x,t,u=EulerJawny(dt,tmax,dx,xmax,0.1)
	T,X=np.meshgrid(t,x)
	my_col = cm.jet(u/np.amax(u))
	cmap = plt.get_cmap('hot')
	fig=plt.figure()
	ax1=fig.add_subplot(1,2,1,projection='3d')
	ax1.set_xlabel('x')
	ax1.set_ylabel('t')
	ax1.set_zlabel('U(x,t)')
	ax1.plot_surface(X,T,u,cmap=cmap)
	ax1.set_title("wersja 3D")
	ax2=fig.add_subplot(1,2,2)
	ax2.pcolormesh(X,T,u,cmap=cmap)
	ax2.set_title("wersja 2D")
	ax2.set_xlabel('x')
	ax2.set_ylabel('t')
	ax2.plot()
	plt.subplots_adjust(wspace = 0.7,hspace = 10 )
	plt.show()
else:
	print("Warunek stabilności nie spełniony")
