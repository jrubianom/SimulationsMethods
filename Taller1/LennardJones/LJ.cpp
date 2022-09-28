// Simular el movimiento de N granos con gravedad, fricción y choques inelásticos
#include <fstream>
#include <iostream>
#include "Particulas.h"

using namespace std;

void Evolucion(Parametro P);

const double g=0, eps=1.0, r0 = 10, R = 2.5, K = 1e4;
const double Lx=60, Ly=120;
const int Nx=5, Ny=5, N=Nx*Ny;

//-----------  Programa Principal --------------  
int main(void){
  double m0=1, R0=R, kT=0.1, V0=sqrt(2*kT/m0);
  double tmax= 400;
  double tcuadro=tmax/(5*100),dt=1e-3;
  Parametro P;
  P.Init(N,Nx,Ny,Lx,Ly,m0,R0,kT,tmax,dt,g,eps,r0,K,tcuadro,true);
  Evolucion(P);
  return 0;
}

double Energia(Cuerpo *Granos,Parametro P);
double yprom(Cuerpo *Granos,Parametro P);
  
void Evolucion(Parametro P){
  int N = P.N;
  Cuerpo Grano[N+4];
  Colisionador Hertz;
  Hertz.Init(N,P.g,P.eps,P.r0,P.K);
  Crandom ran64(1);
  double m0=P.m0, R0=P.R0, kT=P.kT, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,tmax=P.tmax,tcuadro=P.tcuadro,dt=P.dt;
  double Lx = P.Lx, Ly = P.Ly;
  int Nx = P.Nx, Ny = P.Ny;
  double dx=Lx/(Nx+1), dy=Ly/(2*(Ny+1));
  double Theta;
  bool gif = P.gif;

  ofstream fileanim("2DGas.txt");
  Animacion gasAnim;
  gasAnim.Init(-10,Lx+10,-10,Ly+10,Lx,Ly);

  if(gif) gasAnim.InicieAnimacion(fileanim); //Dibujar

  //Inicializar las paredes
  double Rpared=1e7*Lx, Mpared=1e7*m0;
  //---------------(  x0,       y0,Vx0,Vy0,    m0,    R0, theta0,omega0)
  Grano[N+0].Inicie(Lx/2,Ly+Rpared,  0,  0,Mpared,Rpared); //Pared de arriba
  Grano[N+1].Inicie(Lx/2,  -Rpared,  0,  0,Mpared,Rpared); //Pared de abajo
  Grano[N+2].Inicie(Lx+Rpared,Ly/2,  0,  0,Mpared,Rpared); //Pared derecha
  Grano[N+3].Inicie(  -Rpared,Ly/2,  0,  0,Mpared,Rpared); //Pared izquierda
  //Inicializar las moléculas

  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //--------------------(   x0,   y0,          Vx0,          Vy0, m0,R0,theta0,omega0)
      Grano[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy,V0*cos(Theta),V0*sin(Theta), m0,R0);//OJO
    }

    //Grano[0].Inicie(5,10,0,0, m0*1000,R0);//OJO
    //Grano[1].Inicie(15,10,V0,0, m0,R0);//OJO

  ofstream fileYprom("data.txt");
  ofstream fileVelocidades("velx.txt");
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(gif && tdibujo>tcuadro){
      gasAnim.InicieCuadro(fileanim);
      for(i=0;i<N;i++) Grano[i].Dibujese(fileanim);
      gasAnim.TermineCuadro(fileanim);
      tdibujo=0;
    }

    fileYprom << t << "\t" << Energia(Grano,P)<< "\t" << yprom(Grano,P) << "\t" <<endl;
    //--- Muevase por PEFRL ---
    Hertz.Pefrl_unpaso(Grano,dt);

  }

  fileYprom.close();
  fileVelocidades.close();
}

double Energia(Cuerpo *Granos,Parametro P){
  double sum = 0;
  for(int i = 0; i < P.N; i++){
    sum += Granos[i].Energia(P);
  }
  return sum;
}

double yprom(Cuerpo *Granos,Parametro P){
  double sum = 0;
  for(int i=0; i < P.N; i++){
    sum+= Granos[i].Gety();
  }
  return sum/N;
}
