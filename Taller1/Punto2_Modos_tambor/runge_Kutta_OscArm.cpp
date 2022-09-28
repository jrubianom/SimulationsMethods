#include <iostream>
#include <cmath>

using namespace std;

const double omega=3.0;

double ds_dt(double t, double s, double i){
    return -omega*omega * i;
}


double di_dt(double t, double s, double i){
    return s;
}

void UnPasoDeRungeKutta4(double & t, double & s, double & i, double dt){
    double ds1, ds2, ds3, ds4, di1, di2, di3, di4;
    ds1 = dt* ds_dt(t,s,i);                       di1 = dt* di_dt(t,s,i);
    ds2 = dt* ds_dt(t+dt/2,s+ds1/2,i+di1/2);    di2 = dt* di_dt(t+dt/2,s+ds1/2,i+di1/2);
    ds3 = dt* ds_dt(t+dt/2,s+ds2/2, i+di2/2);   di3 = dt* di_dt(t+dt/2,s+ds2/2, i+di2/2);
    ds4 = dt* ds_dt(t+dt/2,s+ds3, i+di3);       di4 = dt* di_dt(t+dt/2,s+ds3, i+di3);

    s+= (ds1 + 2*ds2 + 2*ds3 + ds4)/6;     i+= (di1 + 2*di2 + 2*di3 + di4)/6;
    t+=dt;
}

int main(){
    double t, x, Dx; double dt=0.01;

    for(t=0,x=1, Dx=1; t<2+dt/2; ){
        cout<<t<<" "<<x<<endl;
        UnPasoDeRungeKutta4(t,x, Dx,dt);
    }
    return 0;
}