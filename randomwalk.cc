#include<iostream>
#include<cstdlib>
#include<time.h>
#include<cmath>

int main()
{
  using namespace std;
  int iseed=time(NULL);
  srand48(iseed);

  cerr<<"Quanti passi? ";
  int n;
  cin>>n;
  double th=0,dx=0,dy=0,sumx=0,sumy=0,r=0;

  for(double i=0;i<n;i++)
    {
      th=drand48()*2*M_PI;     
      dx=cos(th);
      dy=sin(th);
      sumx +=dx;
      sumy +=dy;

      // cout<<sumx<<"  "<<sumy<<endl;//
      cout<<sqrt(i)<<" "<<sqrt(sumx*sumx+sumy*sumy)<<endl;
    }
  return 0;
}
