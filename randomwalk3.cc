#include<iostream>
#include<cstdlib>
#include<time.h>
#include<cmath>

using namespace std;

int main ()
{
  int iseed=time(NULL);
  srand48(iseed);
  int n;
  cerr<<"numero di passi: ";
  cin>>n;

  double sumx=0,sumy=0;

  for (int i=0 ;i<=n;i++)
    {
    
      double dx=2*M_SQRT2*(drand48()-0.5);
      double dy=2*M_SQRT2*(drand48()-0.5);

      sumx +=dx ;
      sumy +=dy;

      cout <<sumx<< "  " <<sumy<< endl;
    }

  return 0;
}

