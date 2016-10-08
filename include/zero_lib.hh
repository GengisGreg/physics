#ifndef ZERO_LIB
#define ZERO_LIB

using namespace std;

void zpmbrac(double (*func)(double x),double &x1,double &x2,int &ifail)
{
  int ntry=100;
  ifail=1;
  double hstep=(x2-x1)/ntry;
  double f1=func(x1);
  x2=x1+hstep;
  double f2=func(x2);
  int ncount=0;
  while((f1*f2>0) && (ncount<ntry))
    {
      x1=x2;
      f1=f2;
      x2=x1+hstep;
      f2=func(x2);
      ncount=ncount+1;
    }
  if(ncount<ntry) ifail=0;
}

void zbrac(double(*func)(double x),double &x1,double &x2)
{
  double gamma=1.6;
  int ntry=50;
  double f1=func(x1);
  double f2=func(x2);
  int ncount=0;
  cout<<"fi & f2 ="<<f1<<" "<<f2<<endl;
  while((f1*f2<0) && (ncount<ntry))
    {cout<<"while"<<endl;
      ncount++;
      if(abs(f1)>abs(f2))
	{cout<<"if"<<endl;
	  x1=x1-gamma*(x2-x1);
	  f1=func(x1);
	}
      else
	{cout<<"else"<<endl;
	  x2=x2+gamma*(x2-x1);
	  f2=func(x2);
	}
    }
  x1;x2;cout<<"end"<<endl;
}

void zbisec(double (*func)(double x),double &x1,double &x2,double &zero,double xacc,double facc)
{
  int ntry=100,ncount=1;
  double dx;	 
  if(func(x1)*func(x2)<0.)
    {
      zero=x1;
      dx=0.5*fabs(x2-x1);
      while((fabs(dx)>xacc) && (fabs(func(zero))>facc) && (ncount<ntry))
	{
	  if(func(x1+dx)*func(x2)<0.)
	    {
	      x1=x1+dx;
	      zero=x1;	     
	    }
	  else 
	    {
	      x2=x2-dx;
	      zero=x2;	    	      
	    }	
	  dx=0.5*dx;
	  ncount++;
	}
    }
  //else cout<<"intervallo errato"<<endl;    
  zero;
}

void znewton(double(*func)(double x),double(*dfunc)(double x),double &zero,double xacc)
{
  int ntry=100,ncount=0;
  do
    {
      ncount++;  
      zero=zero-(func(zero)/dfunc(zero));
    }
  while(fabs(func(zero))>xacc && ncount<ntry);  
  zero;
}

void zfalsi(double(*func)(double x),double &x1,double &x2,double &zero,double xacc)
{
  int ntry=100,ncount=0;
  double x=x1;
  do
    {
      ncount++;
      zero=x;
      x=(x1*func(x2)-x2*func(x1))/(func(x2)-func(x1));
      if(func(x1)*func(x)>0) x1=x;
      else x2=x;
    }
  while(fabs(func(x))>xacc && ncount<ntry);
  //if(ncount<=1) cerr<<"error zfalsi "<<ncount<<" iterazione/i"<<endl; 
 zero;
}





double func(double x)
{
  return (exp(-x)-x);
}
double dfunc(double x)
{
  return (-(exp(-x))-1);
}

#endif
