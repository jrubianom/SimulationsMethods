#include <iostream>
#include <cmath>

using namespace std;

//-----Clase Bessel
class Bessel{
private:
    double lambda;
    double r,dr;
    double R, V;
public:
    void Inicie(double lambda0, double r0, double dr, double R0, double V0);
    double Get_r(){return r;}
    double Get_R(){return R;}
    double dR_dr(double r0, double R0, double V0){return V0;}
    double dV_dr(double r0, double R0, double V0){return -(lambda*lambda*R0 + V0*pow(r0,-1));}
    void UnPasoDeRungeKutta4(void);
};

void Bessel::Inicie(double lambda0, double r0, double dr0, double R0, double V0){
        lambda = lambda0;
        r=r0; dr=dr0; R=R0; V=V0;
}

void Bessel::UnPasoDeRungeKutta4(){
    double dR1, dR2, dR3, dR4,   dV1, dV2, dV3, dV4;
    dR1 = dr* dR_dr(r, R, V);                  dV1 = dr* dV_dr(r, R, V);
    dR2 = dr* dR_dr(r+dr/2, R+dR1/2, V+dV1/2); dV2 = dr* dV_dr(r+dr/2, R+dR1/2, V+dV1/2);
    dR3 = dr* dR_dr(r+dr/2, R+dR2/2, V+dV2/2); dV3 = dr* dV_dr(r+dr/2, R+dR2/2, V+dV2/2);
    dR4 = dr* dR_dr(r+dr, R+dR3, V+dV3);       dV4 = dr* dV_dr(r+dr, R+dR3, V+dV3);
   
    R+= (dR1 + 2*(dR2 + dR3) + dR4)/6;     V+= (dV1 + 2*(dV2 + dV3) + dV4)/6;
    r+=dr;
}


int main(){
    Bessel Tambor;
    double lambda; double r0, R0, dR_dr_0, R; double r, rmax, dr;

    //Condiciones iniciales
    r0=0.01;
    R0=1.0; dR_dr_0=0.0;

    //Parametros
    lambda = 1.0;

    //Parametros de grafica   
     rmax=10.0; dr=0.01;

    //Columnas de datos
    Tambor.Inicie( lambda, r0, dr, R0, dR_dr_0);
    r=Tambor.Get_r();

    for (; r<=rmax; ){
        r=Tambor.Get_r(); R=Tambor.Get_R();
        cout<<r<<"\t"<<R<<endl;
        Tambor.UnPasoDeRungeKutta4();
    }

    return 0;
}
