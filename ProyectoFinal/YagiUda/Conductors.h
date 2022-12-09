#include<iostream>
#include<cmath>
#include "Vector.h"

using namespace std;

class Parasitic;
class Sources;


class Parasitic{
    private:
        //------//Coordinates
        double x,y,z;
        //------//Dimensions
        double lx,ly,lz; //width,length,high
        //------------
        double SigmaP=0.0125; //Conductivity
        //-----------------
        double beta; //How fast Sigma goes to zero at the interface Conductor-empty space

    public:
        Parasitic(double x0,double y0,double z0,double lx0,double ly0,double lz0,double SigmaP0,double beta0=2);
        friend class LatticeBoltzmann;
        friend class Theo;
        double Sigma(int ix,int iy,int iz);
};

Parasitic::Parasitic(double x0,double y0,double z0,double lx0,double ly0,double lz0,double SigmaP0,double beta0){
    x = x0; y=y0; z=z0;
    lx = lx0; ly=ly0; lz=lz0;
    SigmaP=SigmaP0;
    beta=beta0;
}

double Parasitic::Sigma(int ix,int iy,int iz){
    double ax = x - lx/2, bx = x + lx/2;
    double ay = y - ly/2, by = y + ly/2;
    double az = z - lz/2, bz = z + lz/2;
    double sx = (1-tanh(beta*(ax-ix)))*(1-tanh(beta*(ix-bx)))/4.0;
    double sy = (1-tanh(beta*(ay-iy)))*(1-tanh(beta*(iy-by)))/4.0;
    double sz = (1-tanh(beta*(az-iz)))*(1-tanh(beta*(iz-bz)))/4.0;
    return SigmaP*sx*sy*sz;
}



class Sources{
    private:
        //------//Coordinates
        int x,y,z;
        //------//Dimensions
        double lx,ly,lz; //width,length,high
        //------------
        double JS=0.0001; //Current
        //-----------------
        double alphaS,beta; //How fast Sigma goes to zero at the interface Conductor-empty space
        double omegaS,lambdaS,CS,kS;
        bool dimension3D;

    public:
        Sources(int x0,int y0,int z0,double lx0,double ly0,double lz0,
                double lambdaS0,double omegaS0,double JS0,
                double beta0,double alpha0=2,bool dimension3D0=false);
        friend class LatticeBoltzmann;
        friend class Theo;
        vector3D JSource(int ix,int iy,int iz,int t);
        double smooth(double a,double b,double z);
};

Sources::Sources(int x0,int y0,int z0,double lx0,double ly0,double lz0,
                 double lambdaS0,double CS0,double JS0,
                 double beta0,double alpha0,bool dimension3D0){
    x = x0; y=y0; z=z0;
    lx = lx0; ly=ly0; lz=lz0;
    lambdaS=lambdaS0;CS=CS0;omegaS=2*M_PI*CS/lambdaS;kS = omegaS/CS;
    JS = JS0;
    alphaS=alpha0; beta=beta0;
    dimension3D=dimension3D0;

}


vector3D Sources::JSource(int ix,int iy,int iz,int t){
    vector3D Jprima0;
    double Jz,J0p;
    if(dimension3D){
        J0p = JS*smooth(-lz/2,lz/2,iz-z);
        Jz=J0p*exp(-alphaS*(pow((ix-x),2)+pow((iy-y),2)))*sin(omegaS*t)*cos(kS*abs(iz-z));
    }
    else
        Jz=JS*exp(-alphaS*( pow((ix-x),2) + pow((iy-y),2))) *sin(omegaS*t);

    Jprima0.cargue(0,0,Jz);
    return Jprima0;
}

double Sources::smooth(double a,double b, double z){
  return tanh(beta*(z-a)) - tanh(beta*(z-b));
}
