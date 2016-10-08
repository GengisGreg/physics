#include<iostream>
#include<cstdlib>
#include<time.h>
#include<cmath>
#include<math.h>



int main ()
{
using namespace std;

  int iseed=time(NULL);
  srand48(iseed);

  double x,y;
  double R_i , R_medio;
  double sum=0;

  for(int N=100 ; N<= 5000 ; N=N+100){
    for(int j=1 ; j<=100 ; j++){ 
      x=0;
      y=0;
      for (int i = 0 ; i <=N ; i++)
	{
	
	double dx =2*M_SQRT2*(drand48()-0.5);
	double dy =2*M_SQRT2*(drand48()-0.5);
	
	x +=dx;
	y +=dy;
	
      }
      
      R_i=sqrt(x*x+y*y);
      sum +=R_i;
      
    }
    
    R_medio=(1/(double)N)*sum;
    double n=sqrt(N);

    cout <<n<<"  "<<R_medio<< endl;
    
  }

  return 0;
}

