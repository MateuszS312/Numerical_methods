
#include <stdint.h>
#include <iostream>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
double calculate_H(std:: vector<float> &psi,std::vector<float> &x,  float omega,  float kappa, float dx, int kk,float t)
{
	//std::cout<<(psi[kk+1]+psi[kk-1]-2*psi[kk])/(dx*dx)<<std::endl;
	return -0.5*(psi[kk+1]+psi[kk-1]-2*psi[kk])/(dx*dx) + kappa*(x[kk]-0.5)*psi[kk]*std::sin(omega * t);
}
float calculateN(int N, std::vector<float> &psi_I, std::vector<float> &psi_R)
{
	float sum=0;
	for(int ii=0;ii<N+1;ii++)
	{
		sum+=(psi_I[ii]*psi_I[ii]+psi_R[ii]*psi_R[ii]);
	}
	return sum/N;
}
float calculateX(int N, std::vector<float> &psi_I, std::vector<float> &psi_R,std::vector<float> &xk)
{
	float sum=0;
	for(int ii=0;ii<N+1;ii++)
	{
		sum+=xk[ii]*(psi_I[ii]*psi_I[ii]+psi_R[ii]*psi_R[ii]);
	}
	return sum/N;
}
double calculateE(int N, std::vector<float> psi_I, std::vector<float> psi_R,std::vector<float> H_r,std::vector<float> H_i) 
{
	
    double sum = 0;
    for(int k = 0; k < N+1; k++) {
        sum += psi_R[k]*H_r[k] + psi_I[k]*H_i[k];
    }
    return sum/N;
}



int main(int argc,char *arg[])
{
	if(argc==1)
	{
		std::cout<<"brak liczby"<<std::endl;
		return 1;
	}
	float multiplier=atoi(arg[1])/400.0;
	const int N=100;
	const int n=1;
	//double L = 1;
	//	double K_0 = 1;
	//double omega = 3 * M_PI * M_PI / 2;
	double omega = multiplier*4.0 * M_PI * M_PI / 2.0;
	double kappa = 6;
	double dtau = 0.0001;
	double steps = 500000;
	double dx=1.0/N;
	//std::cout<<dx<<std::endl;
	std::vector<float> xk(N+1);
	for(int ii=0;ii<N+1;ii++)
	{
		xk[ii]=(((float)ii)/N);
	}
	//for(auto val:xk)
	//	std::cout<<"xk: "<<val<<std::endl;
	std::vector<float> psi_R(N+1);
	std::vector<float> psi_I(N+1);
	for(int ii=0;ii<N+1;ii++)
	{
		psi_R[ii]=(std::sqrt(2) * std::sin(n * M_PI * xk[ii]));
		psi_I[ii]=(0);
	}
	/*for(auto val:psi_R)
	{
		std::cout<<val<<std::endl;
	}
	*/
	/*for(auto val:psi_I)
	{
		std::cout<<val<<std::endl;
	}*/
	
	std::vector<float> H_R(N+1);
	std::vector<float> H_I(N+1);
	H_R[0]=0;
	H_R[N]=0;
	H_I[0]=0;
	H_I[N]=0;
	for( int ii=1;ii<N;ii++)
	{
		H_R[ii]=calculate_H( psi_R, xk, omega,kappa,dx,ii,0);
		H_I[ii]=calculate_H( psi_I, xk, omega,kappa,dx,ii,0);
	}
	/*for(auto val:H_R)
	{
		std::cout<<val<<std::endl;
	}*/
	std::string n1="out_norm";
	std::string n2="out_x";
	std::string n3="out_e";
	std::ofstream out_norm(((n1.append(arg[1])).append(".dat")));
    std::ofstream out_x(((n2.append(arg[1])).append(".dat")));
	std::ofstream out_e(((n3.append(arg[1])).append(".dat")));
	double tau=0;
	for(int ii=0;ii<steps;ii++)
	{
		for(int kk=0;kk<N+1;kk++)
		{
			psi_R[kk]+=(dtau/2)*H_I[kk];
			//std::cout<<kk<<": "<<psi_R[kk]<<std::endl;
		}
		
		tau+=dtau/2;
		for(int kk=1;kk<N;kk++)
		{
			H_R[kk]=calculate_H( psi_R, xk, omega,kappa,dx,kk,tau);
		}
	
		for(int kk=0;kk<N+1;kk++)
		{
			psi_I[kk]-=dtau*H_R[kk];
			
		}
	
		tau+=dtau/2;
		for(int kk=1;kk<N;kk++)
		{
			H_I[kk]=calculate_H( psi_I, xk, omega,kappa,dx,kk,tau);
		}
		
		
		for(int kk=0;kk<N+1;kk++)
		{
			psi_R[kk]+=dtau/2*H_I[kk];
		}
		
		 if(ii % 2 == 0) 
		 {
			 for(int kk=1;kk<N;kk++)
			{
				H_R[kk]=calculate_H( psi_R, xk, omega,kappa,dx,kk,tau);
			}
			out_norm<<ii<<"\t"<<calculateN(N,psi_I,psi_R)<<std::endl;
			out_x<<ii<<"\t"<<calculateX(N,psi_I,psi_R,xk)<<std::endl;
			out_e<<ii<<"\t"<<calculateE(N,psi_I,psi_R,H_R,H_I)<<std::endl;

		}
	}
	out_norm.close();
    out_x.close();
    out_e.close();
	
return 0;
}
