
#include <iostream>
#include <cmath>
#include <fstream>
#include "Random64.h"
using namespace std;

const int Nx=256 , Ny=256, Q=4;
const double p=0.25,p0=0.25;

//-------Clase LatticeGas

class LatticeGas{
    private:
        int V[Q]; //V[i] i = 0 (up), i=1 (right),i=2 (down), i=3 (left)
        double f[Nx*Ny][Q],fnew[Nx*Ny][Q];  //n[ix][i]
    public:
        LatticeGas();
        void Borrese();
        void Inicie(int N,double mux,double muy,double sigmax,double sigmay);
        void show();
        int index(int ix,int iy){return ((ix+Nx)%Nx)*Ny+((iy+Ny)%Ny);};
        void Colisione();
        void Adveccione();
        double rho(int ix,int iy);
        double varianza();
        void Plotrho();
};

LatticeGas::LatticeGas(){
    V[0]=1;V[1]=1,V[2]=-1,V[3]=-1;
}


void LatticeGas::Borrese(){
    int ci;
    for(int ix=0;ix<Nx;ix++){
        for(int iy=0;iy < Ny;iy++){
            ci=index(ix,iy);
            for(int i=0;i<Q;i++){
                f[ci][i]=fnew[ci][i]=0;
            }
        }
    }
}

void LatticeGas::Inicie(int N, double mux,double muy, double sigmax,double sigmay){
    int ix,iy,i,ci;
    double rho0;
    for(ix=0;ix < Nx; ix++)
        for(iy=0;iy < Ny; iy++){
            rho0= (1.0/(sigmax*sqrt(2*M_PI))*exp(-0.5*pow((ix-mux)/sigmax,2)))*( 1.0/(sigmay*sqrt(2*M_PI))*exp(-0.5*pow((iy-muy)/sigmay,2)));
            if(rho0>0.0001) cout<<rho0/Q<<endl;
            ci=index(ix,iy);
            for(i=0;i < Q; i++)
                f[ci][i] = rho0/Q;
        }
}



void LatticeGas::Colisione(){
    int ci;// current_index
    for(int ix=0;ix<Nx;ix++){//Para cada celda
        for(int iy=0;iy < Ny; iy++){
            ci = index(ix,iy);
            fnew[ci][0] = p0*f[ci][0] + p*f[ci][1] + p*f[ci][3] + (1-2*p-p0)*f[ci][2];
            fnew[ci][2] = p0*f[ci][2] + p*f[ci][1] + p*f[ci][3] + (1-2*p-p0)*f[ci][0];
            fnew[ci][1] = p0*f[ci][1] + p*f[ci][0] + p*f[ci][2] + (1-2*p-p0)*f[ci][3];
            fnew[ci][3] = p0*f[ci][3] + p*f[ci][0] + p*f[ci][2] + (1-2*p-p0)*f[ci][1];
        }
    }
}


void LatticeGas::Adveccione(){
    int ci,ni; //current_index, new index
    for(int ix=0;ix<Nx;ix++)
        for(int iy=0;iy < Ny; iy++){
            ci=index(ix,iy);
            for(int i =0;i<Q;i++)
                if(i%2 == 0){
                    ni=index(ix,iy+V[i]);
                    f[ni][i]=fnew[ci][i];
                }
                else{
                    ni=index(ix+V[i],iy);
                    f[ni][i]=fnew[ci][i];
                }
        }
}

double LatticeGas::rho(int ix,int iy){
    double sum = 0;
    int ci = index(ix,iy);
    for(int i = 0; i < Q; i++){
        sum+=f[ci][i];
    }
    return sum;
}

double LatticeGas::varianza(){
    int ix,iy; double N,Xprom,Yprom,Sigma2;
    for(N=0,ix=0;ix<Nx;ix++)
        for(iy=0;iy<Ny; iy++){
            N+=rho(ix,iy);
            //cout << rho(ix,iy) << endl;
        }

    for(Xprom=0,Yprom=0,ix=0;ix<Nx;ix++)
        for(iy=0;iy<Ny; iy++){
            Xprom+=ix*rho(ix,iy);
            Yprom+=iy*rho(ix,iy);
        }
    Xprom/=N;
    Yprom/=N;
    for(Sigma2=0,ix=0;ix<Nx;ix++)
         for(iy=0;iy<Ny;iy++)
             Sigma2+=((pow(ix-Xprom,2)+pow(iy-Yprom,2))*rho(ix,iy));
    Sigma2/=N;

    return Sigma2;
}


//------------Programa Principal-------------

int main(){

    ofstream Varianza("VarianzaVsT.txt");
    int N = 2400, tmax = 400;
    double sigma = 16;
    LatticeGas Diffusion;
    Diffusion.Borrese();
    Diffusion.Inicie(N,int(Nx*0.5),int(Ny*0.5),sigma,sigma);


    for(int t=0; t < tmax; t++){
        Diffusion.Colisione();
        Diffusion.Adveccione();
        Varianza << t << " " << Diffusion.varianza() << endl;
    }
    Varianza.close();
    return 0;

}
