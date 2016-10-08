#ifndef RK4
#define RK4

#include<iostream>
#include<iomanip>
#include<math.h>

void rk4 (double &t, double h, 
          void (*ydot) (double t, double *y, double *dy),
          double *y, 
          int neq )
{
  //
  // FISICA NUMERICA (GZP)
  // rk4s.cc - Program to solve a system of first order differential 
  //           equation. Initial value problem. 
  //           4 points Runge-Kupta
  //
  double fk1[neq], fk2[neq], fk3[neq], fk4[neq];
  double  yt[neq],  dy[neq];
  int i;

  double h12=h/2.;
  ydot(t,y,dy);
  for (i=0; i<neq; i++)
   {
     fk1[i]=h*dy[i];
     yt[i]=y[i]+fk1[i]*0.5;
   }
  ydot(t+h12,yt,dy);
  for (i=0; i<neq; i++)
   {
     fk2[i]=h*dy[i] ;
     yt[i]=y[i]+fk2[i]*0.5 ;
   }
  ydot(t+h12,yt,dy);
  for (i=0; i<neq; i++)
   {
     fk3[i]=h*dy[i];
     yt[i]=y[i]+fk3[i];
   }
  ydot(t+h,yt,dy);
  for (i=0; i<neq; i++)
   {   
     fk4[i]=h*dy[i];
     y[i]=y[i]+(fk1[i]+2.*fk2[i]+2.*fk3[i]+fk4[i])/6.;
   }
  t=t+h;
}
#endif
