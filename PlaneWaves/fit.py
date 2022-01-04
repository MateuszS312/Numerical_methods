
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

def func(x,M,Gamma,C1,C2):
		return (1/(2*np.Pi)*Gamma/((x-M)**2+Gamma**2/4)+C1)*C2
		
data = []
with open('e_omega1.txt') as f:
	for line in f.readlines():
		data.append(list(map(float,line.split())))
data=np.array(data)
xdata=data[:,0]
ydata=data[:,1]
bounds=((14.7,0.01,0,1),(14.9,1,8,20))
popt,pcov=curve_fit(func,xdata,ydata,bounds=bounds)
print(popt)
M,Gamma,C1,C2=popt

x=np.linspace(min(xdata),max(xdata),10000)
y=func(x,popt)
plt.plot(xdata,ydata,'ro')
plt.plot(x,y)
plt.show()
