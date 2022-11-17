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
  int n[Lx][Ly][Q], nnew[Lx][Ly][Q]; //n[ix][iy][i] indica ocupación en casilla[ix][iy] y dirección[i]
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N, double mu, double sigma, Crandom & ran64);
  void Show(void);
  void GrafiqueRho(void);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double rho(int ix, int iy); //Inline
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
        n[ix][iy][i] = nnew[ix][iy][i] =0;
    }
  }
}

void LatticeGas::Inicie(int N,double mu,double sigma,Crandom & ran64){
  int ix,iy,i;
  while(N>0){
    ix=(int) ran64.gauss(mu,sigma);//Escoger una celda en x al azar;
    iy=(int) ran64.gauss(mu,sigma);//Escoger una celda en y al azar;
    if(ix<0) ix=0; if(ix>Lx-1) ix=Lx-1; //Corregir en los bordes, si es necesario;
    if(iy<0) iy=0; if(iy>Ly-1) iy=Ly-1; //Corregir en los bordes, si es necesario;
    i=(int) Q*ran64.r(); //Escoger una dirección al azar;
    if(n[ix][iy][i]==0) // Si está vacío
      {n[ix][iy][i]=1; N--;} //pongo una bolita ahí y decrezco N;
  }
}

void LatticeGas::Show(void){
  for (int ix = 0; ix<Lx; ix++){
    for (int iy = 0; iy<Ly; iy++){
      for (int i  = 0; i < Q; i++){
        cout<<n[ix][iy][i];
      }
      cout<<endl;
    }
    cout<<endl;
  }
  cout<<endl;
}

double LatticeGas::rho(int ix, int iy){
  double sum = 0;
  for (int i = 0; i<Q; i++){
    sum+=n[ix][iy][i];
  }
  return sum;
}

void LatticeGas::GrafiqueRho(void){
  for (int ix =0; ix<Lx; ix++)
    for (int iy=0; iy<Ly; iy++)
      cout <<ix<<"\t"<<iy<<"\t"<<rho(ix,iy)<<endl;
}

void LatticeGas::Colisione(Crandom & ran64){
  double out = ran64.r(); 
  for(int ix=0;ix<Lx;ix++){//para cada celda
    for(int iy=0; iy<Ly; iy++){ 
      if(out > p0 && out<= p0+p){//girar a 90º antihorario
        nnew[ix][iy][0]=n[ix][iy][3]; 
	nnew[ix][iy][1]=n[ix][iy][0];
        nnew[ix][iy][2]=n[ix][iy][1];
	nnew[ix][iy][3]=n[ix][iy][2];} //intercambio los contenidos
      else if (out>2*p+p0){//girar a 180º antiorario
        nnew[ix][iy][0]=n[ix][iy][2];
        nnew[ix][iy][1]=n[ix][iy][3];
        nnew[ix][iy][2]=n[ix][iy][0];
        nnew[ix][iy][3]=n[ix][iy][1];} //intercambio los contenidos
      else if (out > p+p0 && out <= 2*p+p0){     //girar 270º antihorario
        nnew[ix][iy][0]=n[ix][iy][1];
        nnew[ix][iy][1]=n[ix][iy][2];
        nnew[ix][iy][2]=n[ix][iy][3];
        nnew[ix][iy][3]=n[ix][iy][0];} //intercambio los contenidos
      else{			      //dejar estático	
	nnew[ix][iy][0]=n[ix][iy][0];
        nnew[ix][iy][1]=n[ix][iy][1];
        nnew[ix][iy][2]=n[ix][iy][2];
        nnew[ix][iy][3]=n[ix][iy][3];} //intercambio los contenidos//no los intercambio
    }
  }
}

void LatticeGas::Adveccione(void){
  for (int ix = 0; ix< Lx; ix++)
    for (int iy = 0; iy<Ly; iy++)
      for (int i  = 0; i < Q; i++)
        n[(ix+Vx[i]+Lx)%Lx][(iy+Vy[i]+Ly)%Ly][i] = nnew[ix][iy][i];
}

double LatticeGas::Varianza(void){
  int ix,iy; double N,Xprom,Yprom,Sigma2;
  //Calcular N
  for(N=0,ix=0; ix<Lx;ix++)
    for (iy=0;iy<Ly; iy++)
      N+=rho(ix,iy);
  //Calcular Xprom Yprom
  for(Xprom=0, Yprom = 0, ix=0 ;ix<Lx;ix++)
    for (iy=0;iy<Ly; iy++){
      Xprom+=ix*rho(ix,iy);
      Yprom+=iy*rho(ix,iy);
    }
  Xprom/=N; Yprom/=N;
  //Calcular Sigma2
  for(Sigma2=0, ix=0 ;ix<Lx;ix++)
    for (iy=0;iy<Ly; iy++){
      Sigma2+=pow(ix-Xprom,2.0)*rho(ix,iy)+pow(iy-Yprom,2.0)*rho(ix,iy);
    }  
  Sigma2/=(N-1);
  return Sigma2;
}
//-------------------------------MAIN--------------------------------------------
int main(void){
  ofstream Varianza("VarianzaVsT.txt");
  LatticeGas Difusion;
  Crandom ran64(1);
  int N=2400; double mu=Lx/2, sigma=16;
  int t, tmax=350;

  Difusion.Borrese();
  Difusion.Inicie(N,mu,sigma,ran64);
  for(t=0;t<tmax;t++){
    Varianza<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
  }
  //Difusion.GrafiqueRho();
  Varianza.close();
  return 0;
}
