
/*
*
* gcc -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas testbessel.c -o testbessel
*/
#include <stdio.h>
#include <gsl/gsl_specfunc.h>
/*
* gsl_sf_bessel_J0()
* gsl_sf_bessel_J1()
* gsl_sf_bessel_Jn()
* gsl_sf_bessel_I0()
* gsl_sf_bessel_I1()
* gsl_sf_bessel_In()
*/
int main()
{
int i,n,nu,MaxNu = 3;
double x,dx,xmin,xmax;
n = 100;
xmin = 0.0; xmax = 10.0;
dx = (xmax - xmin) / n;
for (i = 0; i <= n; i++) {
x = xmin + i * dx;
printf("%f ", x);
for (nu = 0; nu <= MaxNu; nu++)
printf("%f ", gsl_sf_bessel_Jn(nu, x));
printf("\n");
}

 return 0;
} 
