#include<iostream>
#include<iomanip>
#include<cmath>
#include <../include/PHYSICAL_CONST.hh>
#include <../include/zero_lib.hh>

using namespace std;

double V0,a;

double pari(double E)
{
  double alpha=sqrt((E+V0)/Physics::htc2om);
  double beta=sqrt((-E)/Physics::htc2om);
  return alpha*sin(alpha*a)-beta*cos(alpha*a);
}
double derpari(double E)
{
  double alpha=sqrt((E+V0)/Physics::htc2om);
  double beta=sqrt((-E)/Physics::htc2om); 
  return ((1/Physics::htc2om)/alpha)*(sin(alpha*a)+(alpha*a)*cos(alpha*a)+(alpha/beta)*cos(alpha*a)+a*beta*sin(alpha*a));
}
double dispari(double E)
{
  double alpha=sqrt((E+V0)/Physics::htc2om);
  double beta=sqrt((-E)/Physics::htc2om);
  return alpha*cos(alpha*a)+beta*sin(alpha*a);
}
double derdispari(double E)
{
  double alpha=sqrt((E+V0)/Physics::htc2om);
  double beta=sqrt((-E)/Physics::htc2om); 
  return ((1/Physics::htc2om)/alpha)*(cos(alpha*a)-(alpha*a)*sin(alpha*a)-(alpha/beta)*sin(alpha*a)+a*beta*cos(alpha*a));
}

int main()
{
  cout<<setiosflags(ios::fixed);
  cerr<<"inserire V0(>0)(MeV): ";
  cin>>V0;
  cerr<<"inserire semilarghezza buca a : ";
  cin>>a;
  int ifail;
  double x1=-V0,x2=0.;
  int N,k=0,s;
  cerr<<"numero step di energia : ";
  cin>>N;
  double r=V0/N,E=-V0;
  cerr<<"accuratezza in x e y ";
  double xacc,facc;  
  cin>>xacc>>facc;
  double zero=0,av[20],psi[200],step,norm,h,alpha,sum;
  for (int i=1;i<=N;i++)
    {
      E=E+r;
      x1=E-r;
      x2=E;
      if(E>=0.)break;//per il problema del calcolo in prossimita' di zero      
      zpmbrac(pari,x1,x2,ifail);
      if(ifail==0)
	{
	  zbisec(pari,x1,x2,zero,xacc,facc);
	  znewton(pari,derpari,zero,xacc);	  
	  //cerr<<"zero trovato = "<<setw(15)<<setprecision(10)<<zero<<"  "<<"funzione in zero = "<<setw(15)<<setprecision(10)<<pari(zero)<<endl;
	  av[k]=zero;
	  k++;	
	}
    }
  h=a/100;
  for(int l=0;l<=200;l++)
    {
      cout<<-a+l*h;
      for(int j=0;j<k;j++)
	{
	  alpha=sqrt((av[j]+V0)/Physics::htc2om);
	  
	  sum=cos(alpha*(-a))+cos(alpha*a);
	  
	  for(int p=1;p<=200;p++)
	    {
	     if(pow(-1,p)==1)
	       {
		 s=2;
		 sum +=s*cos(alpha*(-a+p*h));
	       }
	     else
	       {
		 s=4;
		 sum +=s*cos(alpha*(-a+p*h));
	       }
	   }
	 norm=sum*h/3;	 
	 step=-a+l*h;	     
	 cout<<"  "<<(1/sqrt(norm*norm))*cos(alpha*step);
	}
      cout<<endl; 
    }
  
  return 0;
}

