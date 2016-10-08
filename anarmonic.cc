#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include "../include/rk4_lib.hh"
#include "../include/PHYSICAL_CONST.hh"

double alpha=2.*pow(10.,4.)*Physics::m_e/pow(Physics::htc,2);
double E;
double pot(double x)
{
  return 0.5*x*x+0.25*x*x*x*x;
}
void ydot(double t, double *y, double *dy)
{
  double phi=y[0];
  double psi=y[1];
  dy[1]=phi;
  dy[0]=alpha*(pot(t)-E)*psi;
}
int main()
{
  using namespace std;
  cerr<<"energia: ";
  cin>>E;
  double y[2];
  y[0]=0.;
  y[1]=0.0001;
  double x=-8.;
  double dx=0.1;
  while(x<8.)
    {
      rk4(x,dx,ydot,y,2);
      cout<<x<<" "<<y[1]<<endl;
      x +=dx;
    } 
  return 0;
}
