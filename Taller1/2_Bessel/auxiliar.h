#include <fstream>
#include <iostream>
#include <cmath>
#include "ZerosBessel.h"
#include <string>
#include "fstream"

class parameters{
    public:
        double beta, gama;
        double lambda;
        double Error = 1e-5;
        double dt,maxTimeRK;
        double t0; 
        //just in case t0 is a singular point

};

double Bisection(double a,double b,parameters p,double (*func)(double ,parameters )){
    double m = (a+b)/2;
    double fm = func(m,p),fa = func(a,p);
    while(b-a >= p.Error){
        if(fa*fm > 0){
            a = m;
            fa = fm;
        }
        else{
            b = m;
        }
        m = (a+b)/2;
        fm = func(m,p);
    }
    return m;
}

void RungeKuttaOneStep(double &x1, double &x2,double &t,
                       const parameters p,
                       double (*f1)(double,double,double,const parameters),
                       double (*f2)(double,double,double,const parameters)){
    double dx11,dx21,dx31,dx41;
    double dx12,dx22,dx32,dx42;
    double dt = p.dt;
    dx11 = dt*f1(t,x1,x2,p); dx12 = dt*f2(t,x1,x2,p);
    dx21 = dt*f1(t+dt/2,x1+dx11/2,x2+dx12/2,p); dx22 = dt*f2(t+dt/2,x1+dx11/2,x2+dx12/2,p);
    dx31 = dt*f1(t+dt/2,x1+dx21/2,x2+dx22/2,p); dx32 = dt*f2(t+dt/2,x1+dx21/2,x2+dx22/2,p);
    dx41 = dt*f1(t+dt,x1+dx31,x2+dx32,p); dx42 = dt*f2(t+dt,x1+dx31,x2+dx32,p);

    x1 += (dx11 + 2*(dx21 + dx31) + dx41)/6;
    x2 += (dx12 + 2*(dx22 + dx32) + dx42)/6;

    t += dt;
}


void RungeKutta(double &x1, double &x2,double &t,
                const parameters p,
                double (*f1)(double,double,double,const parameters),
                double (*f2)(double,double,double,const parameters),
                std::string namefile){
    std::ofstream file(namefile);
    while(t < p.maxTimeRK){
        RungeKuttaOneStep(x1,x2,t,p,f1,f2);
        file << t << "\t" << x1 << "\t" << x2 << "\t" << std::endl;
    }
    file.close();
}

void RungeKutta(double &x1, double &x2,double &t,
                const parameters p,
                double (*f1)(double,double,double,const parameters),
                double (*f2)(double,double,double,const parameters)){
    while(t < p.maxTimeRK){
        RungeKuttaOneStep(x1,x2,t,p,f1,f2);
    }
}
