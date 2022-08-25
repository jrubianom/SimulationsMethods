#include <iostream>
#include <cmath>
#include <sys/types.h>
#include <type_traits>

using namespace std;

class parameters{
    public:
        double beta,gama;
};

//x1 = s, x2 = i
//SIR model :
// s' = - beta*s*i = f1(t,s,i)
// i' = beta*s*i - gamma*i = f2(t,s,i)
// 1 = s + i + r

double f1(double t,double x1,double x2,const parameters p){
    return -p.beta*x1*x2;
}

double f2(double t,double x1,double x2, const parameters p){
    return (p.beta*x1 - p.gama)*x2;
}

void RungeKuttaOneStep(double &x1, double &x2,double &t,
                       double dt,const parameters p){
    double dx11,dx21,dx31,dx41;
    double dx12,dx22,dx32,dx42;

    dx11 = dt*f1(t,x1,x2,p); dx12 = dt*f2(t,x1,x2,p);
    dx21 = dt*f1(t+dt/2,x1+dx11/2,x2+dx12/2,p); dx22 = dt*f2(t+dt/2,x1+dx11/2,x2+dx12/2,p);
    dx31 = dt*f1(t+dt/2,x1+dx21/2,x2+dx22/2,p); dx32 = dt*f2(t+dt/2,x1+dx21/2,x2+dx22/2,p);
    dx41 = dt*f1(t+dt,x1+dx31,x2+dx32,p); dx42 = dt*f2(t+dt,x1+dx31,x2+dx32,p);

    x1 += (dx11 + 2*(dx21 + dx31) + dx41)/6;
    x2 += (dx12 + 2*(dx22 + dx32) + dx42)/6;

    t += dt;
}

int main(){

    double s,i,r,t,dt = 0.01;
    double beta = 0.35, gama = 0.08;
    parameters p;
    p.beta = beta; p.gama = gama;
    //Initial conditions
    s = 0.999; i = 0.001; r = 1 - s -i;
    //Critical time
    double T = 1/gama;
    double tmax = 5*T;
    for(t = 0;t < tmax;){
        RungeKuttaOneStep(s,i,t,dt,p);
        r = 1-(s+i);
        cout << t << "\t" << s << "\t" << i << "\t" << r << endl;
    }
    return 0;
}
