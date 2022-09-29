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


//F(lambda) = R(r = 1,lambda), i.e, 
//gets the value of R at 1 for a given lambda
double F_lambda(double lambda,parameters p){
    p.lambda = lambda;
    double x1 = 0, x2 = 1, r = p.t0;
    RungeKutta(x1,x2,r,p,&f1,&f2);
    return x2;
}

int main(){
    //Recordemos que los parámetros t realmente son r
    parameters p;
    p.t0 = 1e-8; p.dt = 0.1; p.maxTimeRK = 1;
    double dlambda = 0.1;
    for(double lambda = 0.1; lambda < 16; lambda += dlambda){
        cout << lambda << "\t" << 0 << "\t" << F_lambda(lambda,p) << endl;
    }
    return 0;
}
