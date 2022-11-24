#include <iostream>
#include <cmath>
#include <fstream>
#include "Random64.h"


using namespace std;

const int Lx = 256; // Número de Celdas en x
const int Ly = 256; // Número de Celdas en y 
const double p0 = 0.25; //Probabilidad de giro nulo 
const double p  = 0.25; //Probabilidad de giro de 90º


const int Q = 4; //Número de direcciones por casilla



//--------------------- Clase LatticeGas ------------
class LatticeGas{
private:
  int Vx[Q], Vy[Q]; //V[i] i = 0(derecha) i = 1 (arriba) i = 2 (izquierda) i = 3 (abajo)
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; //n[ix][iy][i] indica ocupación en casilla[ix][iy] y dirección[i]
  //i = 0(derecha) i = 1 (arriba) i = 2 (izquierda) i = 3 (abajo)
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N, double muX, double muY, double sigmaX, double sigmaY);
  void Colisione(void);
  void Adveccione(void);
  double rho(int ix, int iy);
  void Show(void);
  void GrafiqueRho(void);
  double Varianza (void);
};

LatticeGas::LatticeGas(void){
  Vx[0] = 1; Vx[1] = 0; Vx[2] = -1; Vx[3] = 0;
  Vy[0] = 0; Vy[1] = 1; Vy[2] = 0; Vy[3] = -1;
}

void LatticeGas::Borrese(void){
  for (int ix=0; ix< Lx; ix++){
    for (int iy = 0; iy < Ly; iy++){
      for (int i = 0; i < Q; i++)
        f[ix][iy][i] = fnew[ix][iy][i] =0;
    }
  }
}

void LatticeGas::Inicie(int N, double muX, double muY, double sigmaX, double sigmaY){
  int ix,iy,i;
  double rho0;
  for(ix=0;ix < Lx; ix++)
    for(iy=0;iy < Ly; iy++){
      rho0= (1.0/(sigmaX*sqrt(2*M_PI))*exp(-0.5*pow((ix-muX)/sigmaX,2)))*( 1.0/(sigmaY*sqrt(2*M_PI))*exp(-0.5*pow((iy-muY)/sigmaY,2)));//2D gaussian
      for (int i  = 0; i < Q; i++){
        f[ix][iy][i]=rho0/Q; //Se distribuye sobre las posibles direcciones
      }
    }
}



	
void LatticeGas::Colisione(void){
  //i = 0(derecha) i = 1 (arriba) i = 2 (izquierda) i = 3 (abajo)
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0; iy<Ly; iy++){ 
      for (int i = 0; i<Q; i++){
        int inew90 = (i+1)%Q, inew180 = (i+2)%Q, inew270 = (i+3)%Q;      
        fnew[ix][iy][i]= p0*f[ix][iy][i]+p*f[ix][iy][inew90]+p*f[ix][iy][inew270]+(1-2*p-p0)*f[ix][iy][inew180];
      }
    }
  }
}

void LatticeGas::Adveccione(void){
  for (int ix = 0; ix< Lx; ix++)
    for (int iy = 0; iy<Ly; iy++)
      for (int i  = 0; i < Q; i++)
        f[(ix+Vx[i]+Lx)%Lx][(iy+Vy[i]+Ly)%Ly][i] = fnew[ix][iy][i];
}




double LatticeGas::rho(int ix, int iy){
  double sum = 0;
  for (int i = 0; i<Q; i++){
    sum+=f[ix][iy][i];
  }
  return sum;
}

void LatticeGas::Show(void){
  for (int ix = 0; ix<Lx; ix++){
    for (int iy = 0; iy<Ly; iy++){
      for (int i  = 0; i < Q; i++){
        cout<<f[ix][iy][i];
      }
      cout<<endl;
    }
    cout<<endl;
  }
  cout<<endl;
}


void LatticeGas::GrafiqueRho(void){
  for (int ix =0; ix<Lx; ix++)
    for (int iy=0; iy<Ly; iy++)
      cout <<ix<<"\t"<<iy<<"\t"<<rho(ix,iy)<<endl;
}


double LatticeGas::Varianza(){
    int ix,iy; double N,Xprom,Yprom,Sigma2;
    for(N=0,ix=0;ix<Lx;ix++)
        for(iy=0;iy<Ly; iy++){
            N+=rho(ix,iy);
        }

    for(Xprom=0,Yprom=0,ix=0;ix<Lx;ix++)
        for(iy=0;iy<Ly; iy++){
            Xprom+=ix*rho(ix,iy);
            Yprom+=iy*rho(ix,iy);
        }
    Xprom/=N;
    Yprom/=N;
    for(Sigma2=0,ix=0;ix<Lx;ix++)
         for(iy=0;iy<Ly;iy++)
             Sigma2+=((pow(ix-Xprom,2)+pow(iy-Yprom,2))*rho(ix,iy));
    Sigma2/=N;

    return Sigma2;
}

//-------------------------------MAIN--------------------------------------------
int main(void){
  ofstream Varianza("VarianzaVsT.txt");
  LatticeGas Difusion;
  int N=2400; 
  int mu = Lx/2;double sigma=16;
  int t, tmax=400;

  Difusion.Borrese();
  Difusion.Inicie(N,mu,mu,sigma,sigma);
  for(t=0;t<tmax;t++){
    Varianza<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
    //Difusion.Show();
  }
  Varianza.close();
  return 0;
}
