#include <iostream>
#include <cmath>
#include "auxiliar.h"

using namespace std;

// The equation r²R'' + rR' + lambda² r² R = 0 is equiavlent to
// the system
// (x1,x2) = (rR',R) with (x1(0),x2(0)) = (0,1)
// x1' = - lambda² r x2 = f1(r,x1,x2)
// x2' = x1/r = f2(r,x1,x2)

double f1(double t,double x1,double x2,const parameters p){
    return -p.lambda*p.lambda*t*x2;
}

double f2(double t,double x1,double x2, const parameters p){
    return x1/t;
}
//F(lambda) = R(r = 1,lambda)
double F_lambda(double lambda,parameters p){
    p.lambda = lambda;
    double x1 = 0, x2 = 1, r = p.t0;
    RungeKutta(x1,x2,r,p,&f1,&f2,false);
    return x2;
}


int main(){
    double zGuess[6] = {2,4,6,10,14,16};
    parameters p;
    p.t0 = 1e-8; p.dt = 0.01; p.maxTimeRK = 1; p.Error = 1e-5;

    double zeroRK, zeroSimpson;
    double a,b; // the usual a,b of the bisection method
    cout << "Zero\t" << "Actual Zero" << "Error" << endl;
    for(int zero_ith = 0; zero_ith < 5; zero_ith++){
        a = zGuess[zero_ith]; b = zGuess[zero_ith+1];
        zeroRK = Bisection(a,b,p,&F_lambda);
        zeroSimpson = ZeroBessel(0,a,b,100);
        cout << zeroRK << "\t" << zeroSimpson <<
            "\t" << abs(zeroRK - zeroSimpson) << endl;
    }
    return 0;
}
