#include<iostream>
#include<iomanip>
#include<math.h>

#include "../include/zero_lib.hh"
#include "../include/rk4_lib.hh"
#include "../include/PHYSICAL_CONST.hh"

using namespace std;

double E,dx,X;
int num_int;
double*psi;

double V(double x) // potenziale
{
  return E-0.5*x*x-0.25*x*x*x*x;
}

void ydot ( double x2,double*y, double*dy)
{
  double costante =pow(10.,4.)*Physics::m_e/(Physics::htc*Physics::htc);
  dy[0]= y[1];
  dy[1]= -2*costante*V(x2)*y[0];
  
}


double D(double Eps)
  
{ 
  E=Eps;
  //matching point
  double tpoint=0;
  double x_1 = 0;
  double d_x = 0.1;
  double x_2=x_1+d_x;
  double f1 = V(x_1);
  double f2 = V(x_2);
  while(f1*f2 > 0 )
    {
      x_1=x_2;      
      x_2=x_1+d_x;
      f1 = f2;
      f2 = V(x_2);
    }
  tpoint = (x_1+x_2)/2;
  
  double psiup[num_int];
  double psidown[num_int];
  
  int m =(int)((((X)-(tpoint))/(dx))+1);
  m=num_int-m;
  
  
  //psiup
  double y[2];
  y[0]=0.;
  y[1]=pow(10.,-6.);
  
  
  double  x=-X;
  for(int i =0;i<num_int;i++ )
    {
      rk4(x,dx,ydot,y,2);
      psiup[i]=y[0];
    }
  
  //psidown
  y[0]=0.;
  y[1]=pow(10.,-6.);
  
  x=X;
  double dh=-dx;
  for(int l=num_int-1; l>0; l--)
    {
      rk4(x,dh,ydot,y,2);
      psidown[l]=y[0];
    }
  
  double rapporto=(psidown[m]/psiup[m]);
  for(int s=0; s<num_int; s++)
    {
      psidown[s]= psidown[s]/rapporto;
    }
  
  //psi totale
  for(int j=0;j<=num_int;j++)
    {
      if ( j<=m)
	psi[j]=psiup[j];
      else 
	psi[j]=psidown[j];
    }
  
  return (psiup[m-1]-psidown[m-1])/fabs(psiup[m]);
  
} // fine D(E)





int main()
{
  
  E=0; 
  X=5;
  dx=0.1;
  num_int=(int)2*X/dx;
  psi=new double[num_int];
  
  double E1=0;
  double dE=0.1;
  int i=0;
  double av[50];
  double ph[50][num_int+1];
  
  while (E1<50){
    
    double E2=E1+dE;
    
    double D1=D(E1);
    double D2=D(E2);
    
    double acc_E=0.000001;
    double acc_D=0.000001;
    
    if(D1*D2<0){
      double zero = zbisec(D, E1, E2, acc_E,acc_D);
      cerr << "\nL'autovalore numero "<< i
	   <<" si trova in E="<< zero <<endl;
      av[i]=zero;
      
      double sum=0.;
      for(int k=0;k<num_int;k++)
	sum=sum+psi[k]*psi[k];
      sum=sum*dx;
      for(int k=1; k<num_int;k++)
	ph[i][k]=psi[k]/sqrt(sum);
      
      E1=E2;      
      i=i+1;}//if
    else
      E1=E2;
    
  }//while
  
  int num_aut=i;
  cerr << " Autovalori" << endl;
  for(int k=0;k<num_aut;k++)
    cerr << k <<"  " << av[k] <<endl;
  
  double x=-X;
  for(int k=0;k<num_int;k++)
    {
      cout << x <<"  ";
      for(int l=0;l < num_aut;l++)
	cout << ph[l][k] << "  ";
      cout << endl;
      x=x+dx;
    }
  
  return 0;
}
