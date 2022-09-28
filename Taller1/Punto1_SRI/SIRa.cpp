#include <iostream>
#include <cmath>

using namespace std;

const double Beta=0.35, Gamma=0.08;

double ds_dt(double t, double s, double i, double r){
    return -Beta *s*i;
}

double di_dt(double t, double s, double i, double r){
    return (Beta*s - Gamma)*i;
}

double dr_dt(double t, double s, double i, double r){
    return Gamma*i;
}

void UnPasoDeRungeKutta4(double & t, double & s, double & i, double & r, double dt){
    double ds1, ds2, ds3, ds4,   di1, di2, di3, di4, dr1, dr2, dr3, dr4;
    ds1 = dt* ds_dt(t,s,i,r);                         di1 = dt* di_dt(t,s,i,r);
    //dr1 = dt* dr_dt(t,s,i,r);
    ds2 = dt* ds_dt(t+dt/2,s+ds1/2,i+di1/2, r+dr1/2); di2 = dt* di_dt(t+dt/2,s+ds1/2,i+di1/2, r+dr1/2);
    //dr2 = dt* dr_dt(t+dt/2,s+ds1/2,i+di1/2, r+dr1/2);
    ds3 = dt* ds_dt(t+dt/2,s+ds2/2, i+di2/2, r+dr2/2); di3 = dt* di_dt(t+dt/2,s+ds2/2, i+di2/2, r+dr2/2);
    //dr3 = dt* dr_dt(t+dt/2,s+ds2/2, i+di2/2, r+dr2/2);
    ds4 = dt* ds_dt(t+dt,s+ds3, i+di3, r+dr3);       di4 = dt* di_dt(t+dt,s+ds3, i+di3, r+dr3);
    //dr4 = dt* dr_dt(t+dt,s+ds3, i+di3, r+dr3);

    s+= (ds1 + 2*(ds2 + ds3) + ds4)/6;     i+= (di1 + 2*(di2 + di3) + di4)/6;
    //r+= (dr1 + 2*(dr2 + dr3) + dr4)/6;
    t+=dt;
}

int main(){
    double t, s, i, r; double dt=0.01;

    //Condiciones iniciales
    t=0.0;
    s=0.999;
    i=0.001;
    r=1.0-(s+i);

    for( ; t<5.0/Gamma;){
        cout<<t<<"\t"<<s<<"\t"<<i<<"\t"<<r<<endl;
        UnPasoDeRungeKutta4(t,s, i, r, dt);
        r=1-(i+s);
    }
    return 0;
}