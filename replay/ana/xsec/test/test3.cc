#include<gsl/gsl_specfunc.h>
#include<iostream>
//#include<cmath>
void test3(){

int i,n,nu,MaxNu = 3;
double x,dx,xmin,xmax;
n = 100;
xmin = 0.0; xmax = 10.0;
dx = (xmax - xmin) / n;
for (i = 0; i <= n; i++) {
x = xmin + i * dx;
//printf("%f ", x);
for (nu = 0; nu <= MaxNu; nu++)

  gsl_sf_bessel_In(nu, x);
 printf("%f ", gsl_sf_bessel_In(nu, x));

 }

// cout<<" cyl_bessel_i "<< cyl_bessel_i(0,0)<<endl;
 

}
