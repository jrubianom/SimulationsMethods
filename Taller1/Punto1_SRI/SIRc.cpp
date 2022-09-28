#include <iostream>
#include <cmath>

using namespace std;

//-----Clase SRI
class SRI{
    private:
    double beta, gamma;
    double t, s, i, r;
    public:
    void Inicie(double beta0, double gamma0, double t0, double s0, double i0, double r0 );
    double ds_dt(double t0, double s0, double i0){return -beta *s0*i0;}
    double di_dt(double t0, double s0, double i0){return (beta*s0 - gamma)*i0;}
    double dr_dt(double t0, double s0, double i0){return gamma*i0;}
    void UnPasoDeRungeKutta4(double dt);
    void Ejecute(double dt, double tinf);
    double Get_s(){return s;}
};

void SRI::Inicie(double beta0, double gamma0, double t0, double s0, double i0, double r0 ){
        beta=beta0; gamma=gamma0;
        t=t0; s=s0; i=i0; r=r0;
}

void SRI::UnPasoDeRungeKutta4( double dt){
    double ds1, ds2, ds3, ds4,   di1, di2, di3, di4, dr1, dr2, dr3, dr4;
    ds1 = dt* ds_dt(t,s,i);                     di1 = dt* di_dt(t,s,i);
    ds2 = dt* ds_dt(t+dt/2,s+ds1/2,i+di1/2);    di2 = dt* di_dt(t+dt/2,s+ds1/2,i+di1/2);
    ds3 = dt* ds_dt(t+dt/2,s+ds2/2, i+di2/2);   di3 = dt* di_dt(t+dt/2,s+ds2/2, i+di2/2);
    ds4 = dt* ds_dt(t+dt,s+ds3, i+di3);         di4 = dt* di_dt(t+dt,s+ds3, i+di3);
   
    s+= (ds1 + 2*(ds2 + ds3) + ds4)/6;          i+= (di1 + 2*(di2 + di3) + di4)/6;
    t+=dt;
}

void SRI::Ejecute( double dt, double tinf){
    for( ; t<tinf;){
        UnPasoDeRungeKutta4(dt);
        r=1-(i+s);
    }
}


int main(){
    double beta, gamma, R0; double t0, s0, i0, r0;
    //Condiciones iniciales para todos
    t0=0.0;
    s0=0.999;
    i0=0.001;
    r0=1.0-(s0+i0);

    //Parametros
    gamma=0.08;
    beta=gamma;
    

    // t infinito y dt
    double T=1/gamma, tinf=1000*T, dt=0.1*T;

    
    double betafinal=6*beta, dbeta=0.01*beta;
    SRI Epidemia;

    //Columnas de datos
    for (beta=0.0; beta<betafinal; beta+=dbeta){
        R0=beta/gamma;
        Epidemia.Inicie( beta, gamma, t0, s0, i0, r0);
        Epidemia.Ejecute(dt, tinf);
        cout<<R0<<"\t"<<Epidemia.Get_s()<<endl;
    }

    return 0;
}
