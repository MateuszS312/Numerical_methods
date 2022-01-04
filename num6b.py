import numpy as np 
import matplotlib.pyplot as plt

def f1(t,x,y,z):
	return 10*(y-x)
	
def f2(t,x,y,z):
	return x*(28-z)-y

def f3(t,x,y,z):
	return x*y-8/3*z

def RungeKutta4(h,x0,y0,z0,tmax): #krok,wartosc poczatkowa,wspolczynniki
	x=x0
	y=y0
	z=z0	
	X=[]
	Y=[]
	Z=[]
	t=np.arange(0,tmax,h)
	for tt in t:
		X.append(x)
		Y.append(y)
		Z.append(z)
		k1=[]
		k2=[]
		k3=[]
		k4=[]
		#k1
		k1.append(f1(t, x, y, z))
		k1.append(f2(t, x, y, z))
		k1.append(f3(t, x, y, z))
		#K2
		k2.append(f1(t+1/2*h, x + 1/2*k1[0]*h, y + 1/2*k1[1]*h, z + 1/2*k1[2]*h))
		k2.append(f2(t+1/2*h, x + 1/2*k1[0]*h, y + 1/2*k1[1]*h, z + 1/2*k1[2]*h))
		k2.append(f3(t+1/2*h, x + 1/2*k1[0]*h, y + 1/2*k1[1]*h, z + 1/2*k1[2]*h))
		#k3
		k3.append(f1(t+1/2*h, x + 1/2*k2[0]*h, y + 1/2*k2[1]*h, z + 1/2*k2[2]*h))
		k3.append(f2(t+1/2*h, x + 1/2*k2[0]*h, y + 1/2*k2[1]*h, z + 1/2*k2[2]*h))
		k3.append(f3(t+1/2*h, x + 1/2*k2[0]*h, y + 1/2*k2[1]*h, z + 1/2*k2[2]*h))
		#k4
		k4.append(f1(t+h, x + k3[0]*h, y + k3[1]*h, z + k3[2]*h))
		k4.append(f2(t+h, x + k3[0]*h, y + k3[1]*h, z + k3[2]*h))
		k4.append(f3(t+h, x + k3[0]*h, y + k3[1]*h, z + k3[2]*h))
		x=x+1/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])*h
		y=y+1/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])*h
		z=z+1/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])*h
	return t,X,Y,Z	
x0=1
y0=1
z0=1
t,X,Y,Z=RungeKutta4(20/500,x0,y0,z0,20)
fig=plt.figure()
ax1=fig.add_subplot(1,2,1,projection='3d')
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Z')
ax1.set_xticks(np.arange(-40,40,10))
ax1.set_yticks(np.arange(-40,40,10))
ax1.set_zticks(np.arange(0,60,10))
ax1.set_title(r"$ t\in<0, 20>$ - 500  jednostek")
ax1.plot(X,Y,zs=Z,zdir='z')

ax2=fig.add_subplot(1,2,2,projection='3d')
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z')
ax2.set_xticks(np.arange(-40,40,10))
ax2.set_yticks(np.arange(-40,40,10))
ax2.set_zticks(np.arange(0,60,10))
ax2.set_title(r"$ t\in<0, 100>$ - 10000 jednostek")
t,X,Y,Z=RungeKutta4(100/10000,x0,y0,z0,100)
ax2.plot(X,Y,zs=Z,zdir='z')
plt.subplots_adjust(wspace = 0.5,hspace = 10 )
plt.show()
