#include <iostream>
#include <cmath>
#include <string>
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
    RungeKutta(x1,x2,r,p,&f1,&f2);
    return x2;
}


int main(){
    double zGuess[6] = {2,4,6,10,14,16};
    double Zeros[5];
    parameters p;
    p.t0 = 1e-8; p.dt = 0.01; p.maxTimeRK = 1; p.Error = 1e-5;

    double zeroRK, zeroSimpson;
    double a,b; // the usual a,b of the bisection method
    
    string modoExacto;
    //Hallar ceros
    cout << "Zero\t" << "Actual Zero" << "Error" << endl;
    for(int zero_ith = 0; zero_ith < 5; zero_ith++){
        a = zGuess[zero_ith]; b = zGuess[zero_ith+1];
        zeroRK = Bisection(a,b,p,&F_lambda);
        zeroSimpson = ZeroBessel(0,a,b,100);
        cout << zeroRK << "\t" << zeroSimpson <<
            "\t" << abs(zeroRK - zeroSimpson) << endl;
        Zeros[zero_ith] = zeroRK;

        modoExacto = "ModoExacto"+to_string(zero_ith)+".txt";
        printBessel(0,0,1,zeroSimpson,p.dt,modoExacto);
    }

    //Guardar Modos normales Funcines
    double x1 = 0, x2 = 1, r = 0.01,dr = 0.001;
    string filename;
    p.maxTimeRK = 1; p.dt = dr;

    for(int zero_ith = 0; zero_ith < 5; zero_ith++){
        p.lambda = Zeros[zero_ith];
        filename = "Modo"+to_string(zero_ith)+".txt";
        RungeKutta(x1,x2,r,p,&f1,&f2,filename);
        x1 = 0; x2 = 1; r = 0.01;
    }
    return 0;
}
