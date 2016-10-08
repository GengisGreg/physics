#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include "../include/zero_lib.hh"
#include "../include/rk4.hh"
#include "../include/PHYSICAL_CONST.hh"

double alpha=2.*pow(10.,4.)*Physics::m_e/pow(Physics::htc,2);
double E,dx,X;
int nint;
double *PSI;
double pot(double x)
{
  return 0.5*x*x+0.25*x*x*x*x-E;
}
void ydot(double t, double *y, double *dy)
{
  dy[1]=y[0];
  dy[0]=alpha*(pot(t))*y[1];
}
double F(double Ee)
{
  E=Ee;
  double tp=0,x1=0,dx=0.1,x2=x1+dx,f1=pot(x1),f2=pot(x2);
  while(f1*f2>0.)
    {
      x1=x2;
      x2=x1+dx;
      f1=f2;
      f2=pot(x2);
    }
  tp=(x1+x2)/2.;

  int m=nint-(int)((X-tp)/dx);
  double psiup[nint],psidown[nint],y[2];
  //trovo psiup//
  y[0]=0.000001;
  y[1]=0.;
  double x=-X;
  for(int i=0;i<nint;i++)
    {
      rk4(x,dx,ydot,y,2);
      psiup[i]=y[1];
    }
  //trovo psidown//
  y[0]=0.000001;
  y[1]=0.;
  x=X;
  for(int j=nint-1;j>0;j--)
    {
      rk4(x,-dx,ydot,y,2);
      psidown[j]=y[1];
    }
  double rapp=(psidown[m]/psiup[m]);
  for(int l=0;l<nint;l++)
     {
       psidown[l]= psidown[l]/rapp;
     }
  //double PSI[nint];
  for(int t=0;t<=nint;t++)
    {
      if (t<=m)
	PSI[t]=psiup[t];
      else 
	PSI[t]=psidown[t];
    }
  return (psiup[m-1]-psidown[m-1])/fabs(psidown[m]);
}
int main()
{
  using namespace std;
  //E=0.;
  double E1=0.;
  double dE=0.1;
  double E2,aval=0.;
  int i=0;
  X=5.;
  dx=0.1;
  nint=(int)(2*X/dx);
  double av[20];
  PSI=new double[nint];
  double psi[20][nint];
  while(E1<50.)
    {
      E2=E1+dE;
      if(F(E1)*F(E2)<0.)
	{
	  zbisec(F,E1,E2,aval,0.0001,0.0001);   
	  av[i]=aval;	
	  double somma=0.;
	  for(int b=0;b<=nint;b++)
	    {
	      somma +=PSI[b]*PSI[b];    
	    }
	  somma=somma*dx;
	  for(int p=1;p<nint;p++)
	    {
	    psi[i][p]=PSI[p]/sqrt(somma);
	    }
	  i++;		 
	  E1=E2;	
	}
      else E1=E2;
    }
  int nav=i;
  cerr << "gli autovalori sono:"<<endl;
  for(int k=0;k<nav;k++)
    cerr <<k<<" "<<av[k]<<endl;
  int w;
  cerr<<"scegli l'autovalore :";
  cin>>w;
  double x=-X;
  for(int r=0;r<nint;r++)
    {
      cout << x <<" "<<psi[w][r]<<" "<<endl;
      x=x+dx;
    }

  cerr<<"calcolo dell'ortogonalita', scegli 2 autovalori :";
  int a,b;
  cin>>a>>b;
  double ris=0.;
  for(int i=1;i<nint;i++)
    {
     ris +=psi[a][i]*psi[b][i];
    }
  ris=ris*dx;
  cerr<<ris<<endl;
  return 0;
}
