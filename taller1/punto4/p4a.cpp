// Simular el movimiento de N granos con gravedad, fricción y choques inelásticos
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double g=9.8, K=2e4;
const int N = 4;
const double L = 0.12, Dmax = 0.03, M = 0.1; 
const double Lx=1.2*L+N*Dmax, Ly=1.2*L;
//------------------------------------------------------
//--------PERFL------------------------------------------
const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- ----------clase cuerpo -----
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau,l; double I;
  int indice;
public:
  void Inicie(double theta0,double omega0,double l0,double m0,double R0,int numero);
  void BorreFuerza(){F.load(0,0,0); tau=0;};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void AdicioneTorque(double tau0){tau+=tau0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double Gettheta(void){return theta;}; //inline
  double Gettau(void){return tau;};//inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double theta0,double omega0,double l0,double m0,double R0,int numero){
  indice = numero;
  theta=theta0; omega=omega0;
  m=m0;  R=R0; l = l0;
  I=(2.0/5.0)*m*R*R;
  double x0 = l*sin(theta)+numero*2*R, y0=-l*cos(theta);
  double vx0 = omega*l*cos(theta), vy0=omega*l*sin(theta);
  r.load(x0,y0,0); V.load(vx0,vy0,0);  
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  theta+=omega*(Coeficiente*dt);
  double x0 = l*sin(theta)+indice*2*R, y0=-l*cos(theta);
  r.load(x0,y0,0.);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  omega+=tau*(Coeficiente*dt/I);
  double vx0 = omega*l*cos(theta), vy0=omega*l*sin(theta);
  V.load(vx0,vy0,0.);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";
}

//----------------- clase Colisionador ------------------
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Grano);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Grano){
  int i,j; vector3D Fg, Tau, aux;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N+4;i++)
    Grano[i].BorreFuerza();
  //--- Sumar el peso ---
  for(i=0;i<N;i++){
    Fg.load(0,-Grano[i].m*g,0);
    Grano[i].AdicioneFuerza(Fg);
  }
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<N;i++)
    for(j=i+1;j<N+4;j++)
      CalculeFuerzaEntre(Grano[i], Grano[j]);
  //---Calcular Torques----------------------------
  for(i=0;i<N;i++){
    aux.load(i*Dmax,0,0);
    Tau = (Grano[i].r-aux)^Grano[i].F;
    Grano[i].AdicioneTorque(Tau.z());
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2){
  vector3D r21=Grano2.r-Grano1.r;
  double d=r21.norm(), s=Grano1.R+Grano2.R-d;
  if(s>0){
    vector3D n=r21*(1.0/d);
    vector3D F2=n*(K*pow(s,1.5));
    Grano2.AdicioneFuerza(F2);   Grano1.AdicioneFuerza(F2*(-1));
  }   
}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'cunanewton.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange["<<-1.2*L<<":"<<Lx<<"]"<<endl;
  cout<<"set yrange["<<-1.3*Ly<<":0]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Grano[N];
  Colisionador Hertz;
  //Crandom ran64(1);
  double m0=100, R0=1.5;  //double kT=10 V0=sqrt(2*kT/m0);
  int i, j,k;
  double omega0=sqrt(g/L), T=2*M_PI/omega0;
  double t,tdibujo,tmax=0.5*T,tcuadro=tmax/1000,dt=1e-3;
  double Theta0=-M_PI/6;
  
  InicieAnimacion(); //Dibujar
  
  //-------molécula 0---------------
  //             (theta0,omega0,l0,m0,R0,numero)
  Grano[0].Inicie(Theta0, 0, L, M, Dmax/2, 0);
  
  //Inicializar las moléculas
  for(j=1;j<N;j++){
     //-------------(theta0,omega0,l0,m0,R0,numero)
     Grano[j].Inicie(0, 0, L, M, Dmax/2, j);
  }
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt/4, tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(k=0;k<N;k++) Grano[k].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }

    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);  

  }   

  
  return 0;
}

  
