#include <iostream>
#include <cmath>
#include "fstream"
#include "/home/live/repos/Files/Random64.h"

using namespace std;

const int Nx=256,Ny=256;

const double p=0.25,p0=0.25;

class LatticeGas{
    private:
        int *n,*nnew;
        int Nx,Ny,Q=4;
    public:
        LatticeGas(int Nx0,int Ny0);
        ~LatticeGas(){delete[] n; delete[] nnew;};
        int index(int ix,int iy,int i){return (((ix+Nx)%Nx)*Ny+((iy+Ny)%Ny))*Q+i;};
        int index_nbh(int ix,int iy, int nbh);
        void Init(int N,double mu1,double mu2,double sigma1,
                      double sigma2,Crandom & ran64);
        void Collision(Crandom &ran64);
        void Streaming();
        double rho(int ix,int iy);
        double variance();
        void Print(const char * NameFile);

};

LatticeGas::LatticeGas(int Nx0,int Ny0){
    n = new int[Nx0*Ny0*Q];
    nnew = new int[Nx0*Ny0*Q];
    for(int i=0;i < Nx0*Ny0*Q;i++){
        n[i] = 0;
        nnew[i] = 0;
    }

    Nx=Nx0, Ny=Ny0;
}

void LatticeGas::Init(int N,double mu1,double mu2,double sigma1,
                      double sigma2,Crandom & ran64){
    int ix,iy,i,current_index;
    while(N > 0){
        ix=(int) ran64.gauss(mu1,sigma1);
        iy=(int) ran64.gauss(mu2,sigma2);
        i = (int) Q*ran64.r();
        current_index = index(ix,iy,i);
        if(n[current_index] == 0){
            n[current_index] = 1;
            nnew[current_index] = 1;
            N--;
        }
    }
}


int LatticeGas::index_nbh(int ix,int iy, int nbh){

    ///////      |xx|0|xx|
    //////       |3 |A|1 |
    /////        |xx|2|xx|

    int inbh;
    int north=0,east=1,south= 2, west=3;
    if(nbh == south)
        inbh = index(ix,iy+1,south); //because the cell 0 will  stream toward south
    else if(nbh == west)
        inbh = index(ix+1,iy,west); //because the cell 1 will  stream toward west
    else if(nbh == north)
        inbh = index(ix,iy-1,north); //because the cell 2 will  stream toward north
    else
        inbh = index(ix-1,iy,east);  //because the cell 3 will  stream toward east
    return inbh;
}


void LatticeGas::Collision(Crandom &ran64){
    //0 north, 1 east, 2 south, 3 west
    double outp; //output probability exp
    int ci0; //current_index at north
    int ci; //current index
    for(int ix=0;ix < Nx;ix++){
        for(int iy=0;iy < Ny; iy++){
            ci0 = index(ix,iy,0);
            outp = ran64.r();
            if(p0 < outp && outp <= p0+p){
                for(int dir=0; dir < Q; dir++){
                    ci = ci0 + (dir+1)%4;
                    nnew[ci0+dir]=n[ci];
                }
            }
            else if(p0+p < outp && outp <= p0+2*p){
                for(int dir=0; dir < Q; dir++){
                    ci = ci0 + (dir+3)%4;
                    nnew[ci0+dir]=n[ci];
                }
            }
            else if((p0+2*p) < outp){
                for(int dir=0; dir < Q; dir++){
                    ci = ci0 + (dir+2)%4;
                    nnew[ci0+dir]=n[ci];
                }
            }
            else
                for(int dir=0; dir < Q; dir++){
                    ci = ci0 + dir;
                    nnew[ci]=n[ci];
                }
        }
    }

}


void LatticeGas::Streaming(){

    int ci; //current_index
    int civ; //index such that n[cell i]_v = nnew[civ]
    for(int ix=0;ix < Nx;ix++){
        for(int iy=0;iy < Ny; iy++){
            for(int i=0; i < Q; i++){
                ci=index(ix,iy,i);
                civ = index_nbh(ix,iy,i);
                n[ci] = nnew[civ];
            }
        }
    }
}

double LatticeGas::rho(int ix,int iy){
    double sum = 0;
    int ci = index(ix,iy,0);
    for(int i = 0; i < Q; i++){
        sum+=n[ci+i];
    }
    return sum;
}

double LatticeGas::variance(){
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

void LatticeGas::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix,iy;
  for(ix=0;ix<Nx;ix++){
    for(iy=0;iy<Ny;iy++){
      rho0=rho(ix,iy);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

int main(){

    int N = 2400, tmax = 400;
    double sigma = 16;
    Crandom ran64(1);
    LatticeGas Diffusion(Nx,Ny);
    Diffusion.Init(N,int(Nx*0.5),int(Ny*0.5),sigma,sigma,ran64);

    for(int t=0; t < tmax; t++){
        Diffusion.Collision(ran64);
        Diffusion.Streaming();
        cout << t << " " << Diffusion.variance() << endl;
    }
    Diffusion.Print("Plot.dat");
    return 0;
}





/*int ix=3,iy=1;
     for(int i = 0; i < 4; i++)
         cout << "|" << Diffusion.index_nbh(ix,iy,i) << "|";



    return 0;
    for(int iy=Ny-1;iy>-1;iy--){
        for(int ix=0;ix<Nx;ix++){
            for(int i = 0; i < 4; i++)
                cout << "|" << Diffusion.index(ix,iy,i) << "|";
            cout << "\t";
        }
        cout << endl;
    }

    return 0;
*/
