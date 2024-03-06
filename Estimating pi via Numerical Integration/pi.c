/*
Ilya Cable
CS 330
Project 3: Estimating pi via Numerical Integration
11/27/2023
*/
/*

To build and run program (in terminal),
navigate to program folder, then enter:

> gcc -o pi pi.c
> ./pi

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#ifndef M_PI 
#define M_PI 3.14159265358979323846264338327950288
#endif

long double trap(long double(*f)(long double), 
		 long double a, long double b, int n) {
  const long double h = (b - a)/n;
  long double result = f(a)+f(b); //get endpoints
  //for trapezoid, every point other than endpoints occurs twice:
  for(long double i = 1L; i*h<b; i++) result += 2L*f(i*h);
  return result*h/2L;//multiplied by h/2
}

long double simpsons(long double(*f)(long double), 
		 long double a, long double b, int n) {
  const long double h = (b - a)/n;
  long double result = f(a)+f(b); //get endpoints
  //In simpson's 1/3, every odd is multiplied by 4 and every even is multiplied by 2:
  //Odd:
  for(long double i = 1L; i*h<b; i+=2L) result+=4L*f(i*h);
  //Even:
  for(long double i = 2L; i*h<b; i+=2L) result+=2L*f(i*h);
  return result*h/3L;//multiplied by h/3
}

long double simpsons38(long double(*f)(long double), 
		 long double a, long double b, int n) {
  const long double h = (b - a)/n;
  long double result = f(a)+f(b); //get endpoints
  //In simpson's 3/8, the middle elements of an interval are both multiplied by 3:
  for(long double i = 1L; i*h<b; i+=3L) result+=3L*(f(i*h)+f((i+1)*h));
  //Overlapping elements:
  for(long double i = 3L; i*h<b; i+=3L) result+=2L*f(i*h);
  return result*h*3L/8L;//multiplied by 3h/8
}

long double booles(long double(*f)(long double), 
		 long double a, long double b, int n) {
  const long double h = (b - a)/n;
  long double result = 7L*(f(a)+f(b)); //get endpoints
  //second and fourth are mult. by 32:
  for(long double i = 1L; i*h<b; i+=2L) result+=32L*f(i*h);
  //Middle element:
  for(long double i = 2L; i*h<b; i+=4L) result+=12L*f(i*h);
  //Overlapping elements:
  for(long double i = 4L; i*h<b; i+=4L) result+=14L*f(i*h);
  return result*h*2L/45L;//multiplied by 2h/45
}

// f(x)
long double f(long double x) {
  return 4.0L/(1.0L + x*x);
}

// absolute value for long double:
long double absld(long double x){
  if(x<0){
    return -x;
  }
  return x;
}

int main(void) {
  printf("N         Trapezoid           Simpson1/3          Simpson3/8          Boole\n");
  //For all values of i, run integration using all methods with n = 12*2^i
  for(int i = 0; i<17; i++){
    int n = 12 * pow(2,i);
    long double error_trap = absld(M_PI-trap(f,0.0,1.0,n));
    long double simp_result = simpsons(f,0.0,1.0,n);
    long double error_simp_thrd = absld(M_PI-simp_result);
    long double error_simp_eight = absld(M_PI-simpsons38(f,0.0,1.0,n));
    long double error_boole = absld(M_PI-booles(f,0.0,1.0,n));

    // My compiler doesn't recognize long double data type for printf.  I cast the long double to double instead.
    //Shouldn't matter since we're only using 10 sig. digits to right of radix point.
    printf("%-10.0d%-20.10Le%-20.10Le%-20.10Le%-20.10Le\n",n,(double)error_trap,(double)error_simp_thrd,(double)error_simp_eight,(double)error_boole);

  }
  return 0;
}