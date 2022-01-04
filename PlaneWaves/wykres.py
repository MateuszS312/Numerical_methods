from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
nums=np.arange(394,404,1)
for val in nums:
	cmd = "./program01"+' '+ str(val) + ' '
	#print('Python script: ', cmd)
	os.system(cmd)
title = r'Wykres $\epsilon$, dla '
title += f'$\omega$='
title += r'$\frac{4\pi^2}{2}$'
dane=[]
data=[]
for val in nums:
	with open("out_e"+str(val)+".dat") as f:
		print("out_e"+str(val)+".dat")
		for line in f.readlines():
			dane.append(list(map(float, line.split())))
		omega = val/400.0*4 * np.pi * np.pi/ 2;
		dane=np.array(dane)
		t=dane[:,0]
		dane=dane[:,1]
		print(dane);
		if val==400:
			plt.figure(figsize=(12,8))
			plt.title(title, fontsize='x-large')
			plt.plot(t, dane)
			plt.grid()
			plt.legend(('E' ), prop={'size': 14})
			plt.xlabel(r'$\tau$', fontsize=16)
			plt.savefig("fig"+f'omega{omega:.2f}'.replace('.','') + '.png')
			plt.show()	
		with open('e_omega.txt','a') as f2:
			data.append([omega,np.max(dane,0)])
		dane=[]
print(data)
def func(x,M,Gamma,C1,C2):
		return (1/(2*np.pi)*Gamma/((x-M)**2+Gamma**2/4)+C1)*C2
data=np.array(data)	
xdata=data[:,0]
ydata=data[:,1]
#bounds=((14.7,0.01,0,1),(14.9,1,8,20))
bounds = ((19.5, 0.01, 0, 1), (19.9, 1, 8, 20))
popt,pcov=curve_fit(func,xdata,ydata,bounds=bounds)
print(popt)
M,Gamma,C1,C2=popt

x=np.linspace(min(xdata),max(xdata),10000)
y=func(x,*popt)#gwiazdka rozpakowuje tablice
plt.title(r'Dopasowanie krzywej Breita-Wignera.$\omega_0 = \frac{4 \pi^2}{2}$', fontsize='x-large')
plt.xlabel(r'$\omega$', fontsize=16)
plt.ylabel(r'$\langle\epsilon\rangle_{max}$', fontsize=16)
plt.plot(xdata, ydata, 'ro')
plt.plot(x, y)
plt.grid()
plt.legend(('Symulacja', 'Krzywa Breita-Wignera'), prop={'size': 14})
plt.savefig("fig"+f'dopasowanie{omega:.2f}'.replace('.','') + '.png')
plt.show()		
