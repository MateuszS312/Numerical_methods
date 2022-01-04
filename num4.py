import numpy as np 
import matplotlib.pyplot as plt

#Macierz
A=np.array([[6,5,-5],[2,6,-2],[2,5,-1]])
#print(A)
x0=np.random.rand(A.shape[0])
#print(x0)
y0=x0/np.linalg.norm(x0)
epsilon=1e-10
def eigen_vector(y):
	eigenval=0
	eigenvec=np.zeros(y.shape)
	for ii in range(0,150,1):
		eigenvec= A@y
		eigenvec=eigenvec/np.linalg.norm(eigenvec)
		eigenval=eigenvec.T @ A @eigenvec	
		if(np.all(np.abs(y-eigenvec))<epsilon):
			plt.plot(ii+1,eigenval,'ro')	
			print("Obliczono w:",ii+1,"Iteracji",sep='\t')
			break
		plt.plot(ii+1,eigenval,'g^')
		y=eigenvec
	return eigenvec,eigenval
vector,val=eigen_vector(y0)
numpy_val,numpy_vector=np.linalg.eig(A)
np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)
print('Macierz A:	',A,sep='\n',end='\n\n')
print('Wyznaczony wektor własny:',vector,sep='\n',end='\n\n')
print('Wyznaczona najw. wartosc wlasna:',val,sep='\n',end='\n\n')
print('Wartosci wlasne i odpowiadające im wektory własne wedlug numpy:',sep='\n')
for ii in range(0,3):
	print('Wektor i wartość własna nr. ',ii+1)
	print('Wartość własna:	',numpy_val[ii],'Wektor własny:	',numpy_vector[ii],sep='\n',end='\n\n')
plt.grid(True)
plt.xlabel("Iteration")
plt.ylabel("Eigen value")
plt.show()


