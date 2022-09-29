#include <fstream>
#include <iostream>
#include <cmath>
#include "auxiliar.h"

using namespace std;

// The equation r²R'' + rR' + lambda² r² R = 0 is equiavlent to
// the system
// (x1,x2) = (rR',R) with (x1(0),x2(0)) = (0,1)
// x1' = - lambda² r x2
// x2' = x1/r

double f1(double t,double x1,double x2,const parameters p){
    return -p.lambda*p.lambda*t*x2;
}

double f2(double t,double x1,double x2, const parameters p){
    return x1/t;
}

int main(){
    // (x1,x2) = (rR',R), r inicial, dr
    double x1 = 0, x2 = 1, r = 0.01,dr = 0.1;
    //Parámetros: lambda, maxTimeRK = rmax, dt=dr
    parameters p;
    p.lambda = 1; p.maxTimeRK = 10; p.dt = dr;
    //RungeKutta(...) escribe directamente sobre Lambda1.txt
    string filename = "Lambda1.txt";
    RungeKutta(x1,x2,r,p,&f1,&f2,filename);
    return 0;
}
