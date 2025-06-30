#include <stdlib.h>
#include <stdio.h>
#include<iostream>
#include<cmath>
#include<ctime>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;

const int Q = 9;
const int NX = 200;
const int NY = 100;

int e[Q][2] = {{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q] = {4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double rho1[NX+1][NY+1],rho2[NX+1][NY+1],u[NX+1][NY+1][2],f1[NX+1][NY+1][Q],f2[NX+1][NY+1][Q],F1[NX+1][NY+1][Q],F2[NX+1][NY+1][Q];
double tau_f1,tau_f2,dx,dy,dt,niu1,G,fai1[NX+1][NY+1],fai2[NX+1][NY+1],F1x[NX+1][NY+1],F1y[NX+1][NY+1],F2x[NX+1][NY+1],F2y[NX+1][NY+1],totalmass0,totalmass;
int i,j,k,ip,jp,m;
int range,minn,r;
double Force[Q],ux,uy,Fx,Fy;
double Utemp1,Utemp2;
double rho0;
double p_in,p_out;
//
double u1eq[NX+1][NY+1][2],u2eq[NX+1][NY+1][2];
double u_[NX+1][NY+1][2],u0_[NX+1][NY+1][2],u1_[NX+1][NY+1][2],u2_[NX+1][NY+1][2];
double rho[NX+1][NY+1],p[NX+1][NY+1];
int R;
double error;
double Gf,G_1,G_2;
double Gs1,Gs2;
int flag[NX+1][NY+1];
double S1x[NX+1][NY+1],S1y[NX+1][NY+1],S2x[NX+1][NY+1],S2y[NX+1][NY+1];


double feq(int k,double rho, double u1eq[2]);
void output(int m);
void Error();

double feq(int k,double rho, double u,double v)
{
	double eu,uv,feq;
	eu = (e[k][0]*u+e[k][1]*v);
	uv = (u*u+v*v);
	feq = w[k]*rho*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
	return feq;
}

void Error()
   {
   	double temp1,temp2;
   	temp1=0;
   	temp2=0;
   	for(i=1;i<NX;i++)
   	   for(j=1;j<NY;j++)
   	   {
   	   	temp1 +=(
		(u_[i][j][0]-u0_[i][j][0])*(u_[i][j][0]-u0_[i][j][0])+(u_[i][j][1]-u0_[i][j][1])*(u_[i][j][1]-u0_[i][j][1]));
		  temp2 +=(u_[i][j][0]*u_[i][j][0]+u_[i][j][1]*u_[i][j][1]);
		  }
		  temp1=sqrt(temp1);
		  temp2=sqrt(temp2);
		    error=temp1/(temp2+1e-30);
   }

void output(int m)
	{
		
		cout<<"iterations = "<<"  "<<m<<" "<<"The max relative error of uv is:"<<setiosflags(ios::scientific)<<error<<endl;
		
		
		ostringstream name;
		if(m==0)
		{
			name<<"Phase_sepa0000"<<m<<".dat";	
			
		}	
		
		if(m>1&&m<100)
		{
			name<<"Phase_sepa000"<<m<<".dat";	
			
		}
		
		if(m>=100&&m<1000)
		{
			name<<"Phase_sepa00"<<m<<".dat";	
			
		}
		
		if(m>=1000&&m<10000)
		{
			name<<"Phase_sepa0"<<m<<".dat";	
			
		}
		
		if(m>=10000&&m<1000000)
		{
			name<<"Phase_sepa"<<m<<".dat";	
			
		}
		
	  ofstream out(name.str().c_str());
	  out<< "Title= \"LBM Phase_sepa\"\n" << "VARIABLES=\"X\",\"Y\",\"rho1\",\"rho\",\"p\",\"rho2\"\n" << "ZONE T=\"BOX\",I=" << NX+1 << ",J=" << NY+1 << ",F=POINT" << endl;
	  for(j=0;j<=NY;j++)
	     for(i=0;i<=NX;i++)
	     {
	     	out<<double(i)<<" "<<double(j)<<" "<<rho1[i][j]<<" "<<rho[i][j]<<" "<<p[i][j]<<" "<<rho2[i][j]<<endl;
		 }
		
		
		
	}



int main()
{
	
	dx = 1.0;
	dy = 1.0;
	dt = dx;
	rho0 = 1.0;
	tau_f1 = 1.0;
	tau_f2 = 1.0;
	G = 0.225;    //Kang(2002)D2Q9(G is nine times larger than our Gc)
	G_1 = -0.15;
	G_2 = G_1; 
	Gf = 7.0;
	R = 25.0;
	Gs1 = 0.3;
	Gs2 = -0.3;    //根据Huang G_ads,k 之间的差异决定了接触角 差异越大，接触角越大 
	
	niu1 = (tau_f1-0.5)/3.0;
	std::cout<<"niu = "<<niu1<<endl;
	
	
	for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			flag[i][j] = 0;	
			
		}	
		
	}	
	
	for(i=0;i<=NX;i++)
	{
		flag[i][0] = 1;
		flag[i][NY] = 1; 
		 
		
	}

	for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			u1eq[i][j][0] = 0.0;
			u1eq[i][j][1] = 0.0;
			u2eq[i][j][0] = 0.0;
			u2eq[i][j][1] = 0.0;
		
		u_[i][j][0] = 0.0;
		u_[i][j][1] = 0.0;
		
		u0_[i][j][0] = 0.0;
		u0_[i][j][1] = 0.0;
		
		F1x[i][j] = 0.0;
		F1y[i][j] = 0.0;
		F2x[i][j] = 0.0;
		F2y[i][j] = 0.0;
		
			
		//	if(hypot(i-double(NX)/2.0,j-double(NY)/2.0)<=R)
		if(hypot(i-double(NX)/2.0,j-R-1)<=R)
		{
		rho1[i][j] = 8.0;
	//	rho2[i][j] = 0.12;
		
		rho2[i][j] = 0.24;    //同样的参数设置下，如果同一个组分下的第二相密度设置很小，就容易扩散 
		}
		
		else
		{
		rho1[i][j] = 0.24;
		
	//	rho1[i][j] = 0.12;
		rho2[i][j] = 8.0;	
			
		}
			
			
		}
		
	}
//	
for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			
			if(flag[i][j]==0)
			{

			fai1[i][j] = rho1[i][j];  //组分1 
			fai2[i][j] = rho2[i][j];  //组分2
			
		//	fai1[i][j] = 1.0-exp(-1.0*rho1[i][j]);  //组分1 
		//	fai2[i][j] = 1.0-exp(-1.0*rho2[i][j]);  //组分2
		
		//	fai1[i][j] = exp(-1.0/rho1[i][j]);  //组分1 
		//	fai2[i][j] = exp(-1.0/rho2[i][j]);  //组分2
			
			}
		} 
			
	}
	
	for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			if(flag[i][j]==0)
			{
			
			F1x[i][j] = 0.0; //组分1力 
			F1y[i][j] = 0.0;
		
			F2x[i][j] = 0.0; //组分2力 
			F2y[i][j] = 0.0; 
			
			
			S1x[i][j] = 0.0; //组分1s力 
			S1y[i][j] = 0.0;
		
			S2x[i][j] = 0.0; //组分2s力 
			S2y[i][j] = 0.0; 
	

			for(k=0;k<Q;k++)
			{
				
			/*	
				if(hypot(e[k][0],e[k][1])>1.1)
				{
					G = Gf/4.0;
				}
				
				else
				{
					G = Gf;
				}
			*/	

				ip = (i+e[k][0]+(NX+1))%(NX+1);
				jp = (j+e[k][1]+(NY+1))%(NY+1);
				
				if(flag[ip][jp]==0)
				{
				 
				F1x[i][j] += -G*fai1[i][j]*w[k]*fai2[ip][jp]*e[k][0];  //组分1力 
				F1y[i][j] += -G*fai1[i][j]*w[k]*fai2[ip][jp]*e[k][1];
			
				F2x[i][j] += -G*fai2[i][j]*w[k]*fai1[ip][jp]*e[k][0];  //组分2力 
				F2y[i][j] += -G*fai2[i][j]*w[k]*fai1[ip][jp]*e[k][1];
				}
				
				
				else
				{
					
				S1x[i][j] += -Gs1*fai1[i][j]*w[k]*e[k][0];  //组分1s力 
				S1y[i][j] += -Gs1*fai1[i][j]*w[k]*e[k][1];
			
				S2x[i][j] += -Gs2*fai2[i][j]*w[k]*e[k][0];  //组分2s力 
				S2y[i][j] += -Gs2*fai2[i][j]*w[k]*e[k][1];	
					
					
					
				}
			//	F1x[i][j] += -G*fai1[i][j]*3.0*w[k]*fai2[ip][jp]*e[k][0];  //组分1力 
			//	F1y[i][j] += -G*fai1[i][j]*3.0*w[k]*fai2[ip][jp]*e[k][1];
			
			//	F2x[i][j] += -G*fai2[i][j]*3.0*w[k]*fai1[ip][jp]*e[k][0];  //组分2力 
			//	F2y[i][j] += -G*fai2[i][j]*3.0*w[k]*fai1[ip][jp]*e[k][1];
				
			//	F1x[i][j] += -G*fai1[i][j]*fai2[ip][jp]*e[k][0];  //组分1力 
			//	F1y[i][j] += -G*fai1[i][j]*fai2[ip][jp]*e[k][1];
			
			//	F2x[i][j] += -G*fai2[i][j]*fai1[ip][jp]*e[k][0];  //组分2力 
			//	F2y[i][j] += -G*fai2[i][j]*fai1[ip][jp]*e[k][1];
			
			
			//	F1x[i][j] += -fai1[i][j]*(G_1*9.0*w[k]*fai1[ip][jp]*e[k][0]+G*9.0*w[k]*fai2[ip][jp]*e[k][0]);  //组分1力 
			//	F1y[i][j] += -fai1[i][j]*(G_1*9.0*w[k]*fai1[ip][jp]*e[k][1]+G*9.0*w[k]*fai2[ip][jp]*e[k][1]);
			
			//	F2x[i][j] += -fai2[i][j]*(G_2*9.0*w[k]*fai2[ip][jp]*e[k][0]+G*9.0*w[k]*fai1[ip][jp]*e[k][0]);  //组分2力 
			//	F2y[i][j] += -fai2[i][j]*(G_2*9.0*w[k]*fai2[ip][jp]*e[k][1]+G*9.0*w[k]*fai1[ip][jp]*e[k][1]);
			
				
			
			
			}
			
	
	}
			
			
		}
		
	}	
	
///	
 

	
//	
	for(i=0;i<=NX;i++)
	{
	 for(j=0;j<=NY;j++)
	 {
	 	for(k=0;k<Q;k++)
	 	{
	
			f1[i][j][k]=feq(k,rho1[i][j],u1eq[i][j][0],u1eq[i][j][1]);
	 		f2[i][j][k]=feq(k,rho2[i][j],u2eq[i][j][0],u2eq[i][j][1]);
	 	
		 }
	  } 
	
	}
	
	
	
	
	

	for(m=0;m<=8e5;m++)
	{
		
//	
for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			if(flag[i][j]==0)
			{

			for(k=0;k<Q;k++)
			{
				//
				//组分1的平衡态速度 
		//u1eq[i][j][0] = u_[i][j][0]+tau_f1*F1x[i][j]/(rho1[i][j]);
		//u1eq[i][j][1] = u_[i][j][1]+tau_f1*F1y[i][j]/(rho1[i][j]);
		
		//组分2的平衡态速度 
		//u2eq[i][j][0] = u_[i][j][0]+tau_f2*F2x[i][j]/(rho2[i][j]);
		//u2eq[i][j][1] = u_[i][j][1]+tau_f2*F2y[i][j]/(rho2[i][j]);
				
				//
				
				
		u1eq[i][j][0] = u_[i][j][0]+tau_f1*(S1x[i][j]+F1x[i][j])/(rho1[i][j]);
		u1eq[i][j][1] = u_[i][j][1]+tau_f1*(S1y[i][j]+F1y[i][j])/(rho1[i][j]);
		
		//组分2的平衡态速度 
		u2eq[i][j][0] = u_[i][j][0]+tau_f2*(S2x[i][j]+F2x[i][j])/(rho2[i][j]);
		u2eq[i][j][1] = u_[i][j][1]+tau_f2*(S2y[i][j]+F2y[i][j])/(rho2[i][j]);
				
				
				
				
				
				F1[i][j][k] = f1[i][j][k]-(f1[i][j][k]-feq(k,rho1[i][j],u1eq[i][j][0],u1eq[i][j][1]))/tau_f1;
				F2[i][j][k] = f2[i][j][k]-(f2[i][j][k]-feq(k,rho2[i][j],u2eq[i][j][0],u2eq[i][j][1]))/tau_f2;
			}
		
			}
		
					
		}	
			
	}	
		
	for(i=0;i<=NX;i++)
	{
		
		
		f1[i][1][2] = F1[i][1][4];
		f1[i][1][5] = F1[i][1][7];
		f1[i][1][6] = F1[i][1][8];
		
		f2[i][1][2] = F2[i][1][4];
		f2[i][1][5] = F2[i][1][7];
		f2[i][1][6] = F2[i][1][8];
		
		
		f1[i][NY-1][4] = F1[i][NY-1][2];
		f1[i][NY-1][7] = F1[i][NY-1][5];
		f1[i][NY-1][8] = F1[i][NY-1][6];
		
		f2[i][NY-1][4] = F2[i][NY-1][2];
		f2[i][NY-1][7] = F2[i][NY-1][5];
		f2[i][NY-1][8] = F2[i][NY-1][6];
		
		
	 }
		
		
		
	for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			if(flag[i][j]==0)
			{

			for(k=0;k<Q;k++)
			{
			
				ip = (i+e[k][0]+(NX+1))%(NX+1);
				jp = (j+e[k][1]+(NY+1))%(NY+1); 
			    f1[ip][jp][k] = F1[i][j][k];
			    f2[ip][jp][k] = F2[i][j][k];
		
			}
			
			}
		}

	}	
	
///
for(i=0;i<=NX;i++)
{
	
	
	for(j=0;j<=NY;j++)
	{
		if(flag[i][j]==0)
		{
		
		
		
		u1_[i][j][0] = 0.0;
		u2_[i][j][0] = 0.0;
		
		u1_[i][j][1] = 0.0;
		u2_[i][j][1] = 0.0;
		
		rho1[i][j] = 0.0;
		rho2[i][j] = 0.0;
		
		
		
		
		for(k=0;k<Q;k++)
		{
			//组分1和组分2 x方向上的速度rhou' 
			u1_[i][j][0] += f1[i][j][k]*e[k][0]; 
			u2_[i][j][0] += f2[i][j][k]*e[k][0];
			
			//组分1和组分2 y方向上的速度rhov' 
			u1_[i][j][1] += f1[i][j][k]*e[k][1]; 
			u2_[i][j][1] += f2[i][j][k]*e[k][1];
			
			//组分1和组分2的密度 
			rho1[i][j] +=f1[i][j][k];
			rho2[i][j] +=f2[i][j][k];
			
		}
		
		rho[i][j] = rho1[i][j]+rho2[i][j];
		
		
		u0_[i][j][0] = u_[i][j][0];
		u0_[i][j][1] = u_[i][j][1]; 
		
		u_[i][j][0] = (u1_[i][j][0]/tau_f1+u2_[i][j][0]/tau_f2)/(rho1[i][j]/tau_f1+rho2[i][j]/tau_f2);
		u_[i][j][1] = (u1_[i][j][1]/tau_f1+u2_[i][j][1]/tau_f2)/(rho1[i][j]/tau_f1+rho2[i][j]/tau_f2);
		
}
				
	}
	

 }


///

Error();
///	
for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			if(flag[i][j]==0)
			{

			fai1[i][j] = rho1[i][j];  //组分1 
			fai2[i][j] = rho2[i][j];  //组分2
			
		//	fai1[i][j] = 1.0-exp(-1.0*rho1[i][j]);  //组分1 
		//	fai2[i][j] = 1.0-exp(-1.0*rho2[i][j]);  //组分2
		
		//	G = 7.0;	
			p[i][j] = rho[i][j]/3.0+G*fai1[i][j]*fai2[i][j]/3.0;
		
		//	G = 27.0;
		//	p[i][j] =  rho[i][j]/3.0+3.0*G*fai1[i][j]*fai2[i][j];
		
		//	G = 3.0;
		//	p[i][j] = rho[i][j]/3.0+3.0*G*fai1[i][j]*fai2[i][j]/2.0;
			
		//	p[i][j] = rho[i][j]/3.0+1.0/6.0*(2.0*G*fai1[i][j]*fai2[i][j]+G_1*pow(fai1[i][j],2)+G_2*pow(fai2[i][j],2));
}
		} 
			
	}
	
	for(i=0;i<=NX;i++)
	{
		for(j=0;j<=NY;j++)
		{
			
			if(flag[i][j]==0)
			{

			F1x[i][j] = 0.0; //组分1力 
			F1y[i][j] = 0.0;
		
			F2x[i][j] = 0.0; //组分2力 
			F2y[i][j] = 0.0; 
		
			S1x[i][j] = 0.0; //组分1s力 
			S1y[i][j] = 0.0;
		
			S2x[i][j] = 0.0; //组分2s力 
			S2y[i][j] = 0.0;
	
			for(k=0;k<Q;k++)
			{
				
				/*
				if(hypot(e[k][0],e[k][1])>1.1)
				{
					G = Gf/4.0;
				}
				
				else
				{
					G = Gf;
				}
				*/
				
				ip = (i+e[k][0]+(NX+1))%(NX+1);
				jp = (j+e[k][1]+(NY+1))%(NY+1);
				
				if(flag[ip][jp]==0)
				{
				 
				F1x[i][j] += -G*fai1[i][j]*w[k]*fai2[ip][jp]*e[k][0];  //组分1力 
				F1y[i][j] += -G*fai1[i][j]*w[k]*fai2[ip][jp]*e[k][1];
			
				F2x[i][j] += -G*fai2[i][j]*w[k]*fai1[ip][jp]*e[k][0];  //组分2力 
				F2y[i][j] += -G*fai2[i][j]*w[k]*fai1[ip][jp]*e[k][1];
				}
				
				else
				{
				
				S1x[i][j] += -Gs1*fai1[i][j]*w[k]*e[k][0];  //组分1s力 
				S1y[i][j] += -Gs1*fai1[i][j]*w[k]*e[k][1];
			
				S2x[i][j] += -Gs2*fai2[i][j]*w[k]*e[k][0];  //组分2s力 
				S2y[i][j] += -Gs2*fai2[i][j]*w[k]*e[k][1];	
						
				}
				
				
			//	F1x[i][j] += -G*fai1[i][j]*3.0*w[k]*fai2[ip][jp]*e[k][0];  //组分1力 
			//	F1y[i][j] += -G*fai1[i][j]*3.0*w[k]*fai2[ip][jp]*e[k][1];
			
			//	F2x[i][j] += -G*fai2[i][j]*3.0*w[k]*fai1[ip][jp]*e[k][0];  //组分2力 
			//	F2y[i][j] += -G*fai2[i][j]*3.0*w[k]*fai1[ip][jp]*e[k][1];
			
			//	F1x[i][j] += -fai1[i][j]*(G_1*9.0*w[k]*fai1[ip][jp]*e[k][0]+G*9.0*w[k]*fai2[ip][jp]*e[k][0]);  //组分1力 
			//	F1y[i][j] += -fai1[i][j]*(G_1*9.0*w[k]*fai1[ip][jp]*e[k][1]+G*9.0*w[k]*fai2[ip][jp]*e[k][1]);
			
			//	F2x[i][j] += -fai2[i][j]*(G_2*9.0*w[k]*fai2[ip][jp]*e[k][0]+G*9.0*w[k]*fai1[ip][jp]*e[k][0]);  //组分2力 
			//	F2y[i][j] += -fai2[i][j]*(G_2*9.0*w[k]*fai2[ip][jp]*e[k][1]+G*9.0*w[k]*fai1[ip][jp]*e[k][1]);
				
			}
			
}
			
			
			
		}
		
	}	
	
///	
 

////	
	
		
		
	
		
		
	
	
			if(m%1000==0)
			{
				
				output(m);
	p_out = p[3*NX/4+R/2][NY/2];
	p_in =  p[NX/2][NY/2];
	std::fstream ll1;
	ll1.open("pin-pout.txt",ios::out|ios::app);	
	ll1<<p_in-p_out<<" "<<3*NX/4+R/2<<" "<<NY/2<<endl;
	ll1.close();
			}
	
	
	
	}
	

	
	return 0;
}
