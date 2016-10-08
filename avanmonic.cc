#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include "../include/zero_lib.hh"
#include "../include/rk4.hh"
#include "../include/PHYSICAL_CONST.hh"

double alpha=2.*pow(10.,4.)*Physics::m_e/pow(Physics::htc,2);
double E,dx,X,psi;
int nint;

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
  return (psiup[m-1]-psidown[m-1])/fabs(psidown[m]);
}
int main()
{
  using namespace std;
  E=0.;
  double dE=0.1;
  X=5.;
  dx=0.1;
  nint=(int)(2*X/dx);
  while(E<50.)
    {
      cout<<E<<" "<<F(E)<<endl;
      E +=dE;
    }
  return 0;
}
