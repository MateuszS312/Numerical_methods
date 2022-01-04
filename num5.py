import numpy as np 
import matplotlib.pyplot as plt
import scipy.linalg as linalg

def aprox(X,Y,m):#x y i stopien wielomianu
	A=np.zeros([m+1,m+1])
	b=[]

	for jj in range(0,2*m+2):
		xji=0
		for val in X:
			xji+=val**jj
		xsum.append(xji)
	for jj in range(0,m+1):
		for ii in range(0,m+1):
			A[ii][jj]=xsum[ii+jj]
		tempb=0
		for x,y in zip(X,Y):
			tempb+=y*x**jj
		b.append(tempb)
	a=np.linalg.solve(A,b)
	print("XSum:",xsum,sep='\n',end='\n\n')
	print("B:",b,sep='\n',end='\n\n')
	print("A:",A,sep='\n',end='\n\n')
	print("a:",a,sep='\n',end='\n\n')
	return a
def poly(coefficients, argument):
	return sum(coef * argument ** power for power, coef in enumerate(coefficients))
#Macierz
X=[1,2,3,4,5]
Y=[11,2,-1,20,-5]
xsum=[]
m=4#stopien wielomianu
if len(X) !=len(Y):
	print("Wiecej wspolrzednych X/Y niz Y/X popraw dane!")
elif len(X)-m<1:
	print("Maksymalny stopień wielomianu musi spelniac relacje ilosc punktów-stopień>=1")
else:
	a=aprox(X,Y,m)
	poly_vec = np.vectorize(poly, excluded=["coefficients"])
	xpoint=np.arange(min(X),max(X),0.001)
	ypoint=poly_vec(argument=xpoint,coefficients=a)
	plt.plot(xpoint,ypoint,'r')
	plt.plot(X,Y,'g^')
	plt.grid(True)
	plt.show()
