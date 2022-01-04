import numpy as np 
import scipy.linalg as linalg

def luDecomp(A):
	
	L, U=np.zeros(A.shape),np.zeros(A.shape)
	suma=0
	N=A.shape[0]
	for ii in range(N):
		L[ii][ii]=1
		for jj in range(ii+1):
			for kk in range(jj):
				suma+=L[jj][kk]*U[kk][ii]
			U[jj][ii]=A[jj][ii]-suma
			#print(U[ii][ii])
			suma=0
		for jj in range(ii+1,N):
			for kk in range(ii):
				suma+=L[jj][kk]*U[kk][ii]
			L[jj][ii]=(A[jj][ii]-suma)/U[ii][ii]
			suma=0	
	return U,L
#definicja macierzy
B=[[9,8,-2,2,-2],[7,-3,-2,7,2],[2,-2,2,-7,6],[4,8,-3,3,-1],[2,2,-1,1,4]]
#B=[[3,2,2],[2,3,2],[2,2,3]]
#B=[[2.0,2.0],[3.0,3.0]]
#B=[[0,0],[0,0]]
A=np.array(B)
if A.shape[0]==A.shape[1]:
	if np.all(A):
		U,L=luDecomp(A)
		np.set_printoptions(precision=2)
		print("Macierze obliczone metoda Doolittle'a:\n")
		print("A: ",A,"U: ", U,"L: ", L, sep='\n\n')
		print("Macierze obliczone przez funkcje wbudowana:")
		p,l,u=linalg.lu(A)
		print("u: ", u,"l: ",l, sep='\n\n')
	else:
		print("Macierz zerowa - zaleznosc spelniona tozsamosciowo")
else:
	print("macierz nie jest kwadratowa\n\n")
	print(A)
