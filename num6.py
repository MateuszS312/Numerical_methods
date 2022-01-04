import numpy as np 
import matplotlib.pyplot as plt

def func(x,y):
	return y-x**2
def analitic(x):
	return 2+2*x+x**2-np.exp(x)
def eulerMethod(h,y0): #krok,wartosc poczatkowa
	X=np.arange(0,3+h,h)
	F=[]
	for val in X:
		F.append(y0)
		y0=y0+h*func(val,y0)
	return F,X
	
def RungeKutta2(h,y0,c1,c2,a2,b21): #krok,wartosc poczatkowa,wspolczynniki
	X=np.arange(0,3+h,h)
	F=[]
	for val in X:
		K1=func(val,y0)
		K2=func(val+a2*h,y0+b21*h*K1)
		F.append(y0)
		y0=y0+(c1*K1+c2*K2)*h
	return F,X	
def RungeKutta4(h,y0,c1,c2,c3,c4,a2,a3,a4,b21,b31,b32,b41,b42,b43): #krok,wartosc poczatkowa,wspolczynniki
	X=np.arange(0,3+h,h)
	F=[]
	for val in X:
		K1=func(val,y0)
		K2=func(val+a2*h,y0+b21*h*K1)
		K3=func(val+a3*h,y0+b31*h*K1+b32*h*K2)
		K4=func(val+a4*h,y0+b41*h*K1+b42*h*K2+b43*h*K3)
		F.append(y0)
		y0=y0+(c1*K1+c2*K2+c3*K3+c4*K4)*h
	return F,X	
	
fig, ax = plt.subplots(1,3)
fig.suptitle("ROZWIAZANIA:", fontsize=16)

x=np.arange(0,3,0.01)
y=[]
for val in x:
	y.append(analitic(val))


for ii in range(0,3):
	ax[ii].set_xlabel("x")
	ax[ii].set_ylabel("y")
	
y0=1
ax[0].set_title('h=0.2')
ax[0].grid(True)
h=0.2
Y,X=eulerMethod(h,y0)
ax[0].plot(X,Y,'--g^')
Y,X=RungeKutta2(h,y0,1/2,1/2,1,1)
ax[0].plot(X,Y,'--ro')
Y,X=RungeKutta4(h,y0,1/6,1/3,1/3,1/6,1/2,1/2,1,1/2,0,1/2,0,0,1)
ax[0].plot(X,Y,'--b+')
ax[0].plot(x,y,'-k')

ax[1].set_title('h=0.02')
ax[1].grid(True)
h=0.2/10
Y,X=eulerMethod(h,y0)
ax[1].plot(X,Y,'--g^')
Y,X=RungeKutta2(h,y0,1/2,1/2,1,1)
ax[1].plot(X,Y,'--ro')
Y,X=RungeKutta4(h,y0,1/6,2/6,2/6,1/6,1/2,1/2,1,1/2,0,1/2,0,0,1)
ax[1].plot(X,Y,'--b+')
ax[1].plot(x,y,'-k')

ax[2].set_title('h=1')
h=5*0.2
Y,X=eulerMethod(h,y0)
ax[2].plot(X,Y,'--g^',label='Euler')
Y,X=RungeKutta2(h,y0,1/2,1/2,1,1)
ax[2].plot(X,Y,'--ro',label='Runge Kutta II')
Y,X=RungeKutta4(h,y0,1/6,1/3,1/3,1/6,1/2,1/2,1,1/2,0,1/2,0,0,1)
ax[2].plot(X,Y,'--b+',label='Runge Kutta IV')
ax[2].plot(x,y,'-k',label='Analityczna')
plt.grid(True)
plt.subplots_adjust(wspace = 0.5,hspace = 10 )
plt.legend()
plt.show()

	
