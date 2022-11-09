#include <iostream>
#include <cmath>
#include "/home/live/repos/Files/Random64.h"
using namespace std;

const int Nx=256,Ny=256,Q=4;
const double p=0.25,p0=0.25;

//-------Clase LatticeGas

class LatticeGas{
    private:
        int V[Q]; //V[i] i = 0 (up), i=1 (right),i=2 (down), i=3 (left)
        int n[Nx*Ny][Q],nnew[Nx*Ny][Q];  //n[ix][i]
    public:
        LatticeGas();
        void Borrese();
        void Inicie(int N,double mux,double muy,double sigmax,double sigmay,Crandom & ran64);
        void show();
        int index(int ix,int iy){return ((ix+Nx)%Nx)*Ny+((iy+Ny)%Ny);};
        void Colisione(Crandom & ran64);
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
                n[ci][i]=nnew[ci][i]=0;
            }
        }
    }
}

void LatticeGas::Inicie(int N, double mux,double muy, double sigmax,double sigmay, Crandom & ran64){
    int ix,iy,i,current_index;
    while(N > 0){
        ix=(int) ran64.gauss(mux,sigmax);
        iy=(int) ran64.gauss(muy,sigmay);
        i = (int) Q*ran64.r();
        current_index = index(ix,iy);
        if(n[current_index][i] == 0){
            n[current_index][i] = 1;
            nnew[current_index][i] = 1;
            N--;
        }
    }
}



void LatticeGas::Colisione(Crandom & ran64){
    double outp;
    int ci;// current_index
    for(int ix=0;ix<Nx;ix++){//Para cada celda
        for(int iy=0;iy < Ny; iy++){
            outp=ran64.r();
            ci = index(ix,iy);
            if(p0 < outp && outp <= p0+p) //Genero numero al azar, sin es mayor que p, volteo
                for(int i=0; i < Q; i++)
                    nnew[ci][i]=n[ci][(i+1)%4];

            else if(p0+p < outp && outp <= p0+2*p)
                for(int i=0; i < Q; i++)
                    nnew[ci][i]=n[ci][(i+3)%4];
            else if((p0+2*p) < outp)
                for(int i=0; i < Q; i++)
                    nnew[ci][i]=n[ci][(i+2)%4];
            else
                for(int i=0; i < Q; i++)
                    nnew[ci][i]=n[ci][i];
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
                    n[ni][i]=nnew[ci][i];
                }
                else{
                    ni=index(ix+V[i],iy);
                    n[ni][i]=nnew[ci][i];
                }
        }
}

double LatticeGas::rho(int ix,int iy){
    double sum = 0;
    int ci = index(ix,iy);
    for(int i = 0; i < Q; i++){
        sum+=n[ci][i];
    }
    return sum;
}

double LatticeGas::varianza(){
    int ix,iy; double N,Xprom,Yprom,Sigma2;
    for(N=0,ix=0;ix<Nx;ix++)
        for(iy=0;iy<Ny; iy++)
            N+=rho(ix,iy);
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
    Sigma2/=(N-1);

    return Sigma2;
}


//------------Programa Principal-------------

int main(){
    int N = 2400, tmax = 350;
    double sigma = 16;
    Crandom ran64(1);
    LatticeGas Diffusion;
    Diffusion.Borrese();
    Diffusion.Inicie(N,int(Nx*0.5),int(Ny*0.5),sigma,sigma,ran64);

    for(int t=0; t < tmax; t++){
        Diffusion.Colisione(ran64);
        Diffusion.Adveccione();
        cout << t << " " << Diffusion.varianza() << endl;
    }
    return 0;

}
