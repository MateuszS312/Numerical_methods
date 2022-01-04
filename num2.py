import numpy as np 
import matplotlib.pyplot as plt

f=lambda x: 3*x**3-4*x**2-1

def bisekcja(aa,bb):
	e=1e-5
	xi=[]
	x0=0
	ii=0
	while abs(f(x0))>e:
		ii+=1
		x0=aa+(bb-aa)/2
		xi.append(x0)
		if (f(x0)*f(aa)>0):
			aa=x0
		elif(f(x0)*f(bb)>0):
			bb=x0
		print("w trakcie obliczen x0= ",x0)
	print("Rozwiazanie x0= ",x0,"\nLiczba kroków: ",ii)
	return x0, xi
a=float(input("Podaj a:"))
print("a= ",a)
b=float(input("Podaj b:"))
print("b= ",b)
if(abs(f(a)*f(b))<1e-5):
	if(abs(f(a))<1e-5):
		print("a to miejsce zerowe")
	elif(abs(f(b))<1e-5):
		print("b to miejsce zerowe")
elif(f(a)*f(b)<0):
	xi=[]
	x0,xi=bisekcja(a,b)
	x=np.linspace(a,b,100)
	x2=np.array(xi)
	y2=f(x2)
	y=f(x)
	plt.plot(x,y)
	plt.plot(x,y,'r',x2,y2,'g^')
	plt.grid(True)
	plt.xlabel("X")
	plt.ylabel("Y")
	plt.show()
else:
	print("a= ",a,"f(a)= ",f(a))
	print("b= ",b,"f(b)= ",f(b))
	print("zly przedzial")
