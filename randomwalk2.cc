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
  double dx=0,dy=0,sumx=0,sumy=0;

  for(int i=0;i<n;i++)
   {
    dx=(2*drand48())-1;
    dy=(2*(rand()%2)-1)*sqrt(1-(dx*dx));
    sumx +=dx;   
    sumy +=dy;
    cout<<sumx<<"  "<<sumy<<endl;
   }
 
return 0;
}
