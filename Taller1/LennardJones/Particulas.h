// Simular el movimiento de N granos con gravedad, fricción y choques inelásticos
#include <iostream>
#include <cmath>
#include "/home/live/repos/Files/vector.h"
#include "Random64.h"
#include "fstream"

using namespace std;

//---- declarar constantes ---
const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;
class Parametro;
class Animacion;

//---- interface e implementacion de clases ----
//--- clase Parametro

class Parametro{
  public:
    //Numero de particulas, y numero de particulas en x, en y para inicializar
    int N,Nx,Ny;
    double Etotal;
    //Dimensiones de la caja
    double Lx,Ly;
    //parametros de las particulas
    double m0,R0,kT;
    //paramestros de las fuerza
    double g;
    double eps, r0; //Lennard-jones
    double K = 1.0e4; //Fuerza paredes
    //tiempos de simulacion
    double tmax,dt;

    //parametros de animacion
    double tcuadro;
    bool gif = false;

    //Inicializar
    void Init(int Nn,int Nnx,int Nny,double Llx,double Lly,double mm0,
              double Rr0,double kTt,double tmaxx,double dtt, double g0,
              double eps0, double r00,double K0,double tcuadro0,bool gif0){
      N = Nn; Nx = Nnx; Ny = Nny; Lx = Llx; Ly=Lly;
      m0=mm0; R0=Rr0; kT=kTt; tmax=tmaxx; dt=dtt;
      g=g0; eps=eps0; r0 = r00;
      tcuadro = tcuadro0; gif=gif0;
      K = K0;
    }
};

//---- clase cuerpo ---
class Cuerpo{
  private:
    vector3D r,V,F; double m,R;
  public:
    void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0);
    void BorreFuerza(){F.load(0,0,0);};
    void AdicioneFuerza(vector3D F0){F+=F0;};
    void Mueva_r(double dt, double Coeficiente);
    void Mueva_V(double dt, double Coeficiente);
    void Dibujese(ofstream &file);
    double Energia(Parametro P);
    double Getx(void){return r.x();}; //inline
    double Gety(void){return r.y();}; //inline
    friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
                    double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
}
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);
}
void Cuerpo::Dibujese(ofstream &file){
  file <<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

double Cuerpo::Energia(Parametro P){
  double E = 0;
  double rr = pow(P.r0/r.norm(),6);
  E += 0.5*m*V.norm2();
  E += P.eps*rr*(rr - 2);
  return E*1e3;
}


//--- clase Colisionador ----
class Colisionador{
  private:
    int N;
    double g=0, eps=1.0, r0 = 10,K;
  public:
    void Init(int N0,double g0,double eps0,double r00,double K0){
      g = g0; eps = eps0; r0=r00;
      N = N0; K = K0;
    }
    void CalculeFuerzas(Cuerpo * Grano);
    void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,int pared);
    void Pefrl_unpaso(Cuerpo *Grano,double dt);
};

void Colisionador::CalculeFuerzas(Cuerpo * Grano){
  int i,j; vector3D Fg;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N+4;i++)
    Grano[i].BorreFuerza();
  //--- Sumar el peso ---
  for(i=0;i<N;i++){
    Fg.load(0,-Grano[i].m*g,0);
    Grano[i].AdicioneFuerza(Fg);
  }
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<N;i++){
    for(j=i+1;j<N+4;j++)
      CalculeFuerzaEntre(Grano[i], Grano[j],j);
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2, int pared){
  vector3D r21=Grano2.r-Grano1.r;
  double d=r21.norm(), s=Grano1.R+Grano2.R-d;
  //si se trata de una particula
  if(pared < N){
    vector3D n=r21*(1.0/d);
    double xx = pow(r0/d,6);
    vector3D F2=n*(12*eps/d)*xx*(xx-1);
    Grano2.AdicioneFuerza(F2);   Grano1.AdicioneFuerza(F2*(-1));
  }
  else if(s>0){ //si se trata de una pared
    vector3D n=r21*(1.0/d);
    vector3D F2=n*(K*pow(s,1.5));
    Grano2.AdicioneFuerza(F2);   Grano1.AdicioneFuerza(F2*(-1));
  }
}

void Colisionador::Pefrl_unpaso(Cuerpo *Grano,double dt){
  //--- Muevase por PEFRL ---
  int i;
  for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
  CalculeFuerzas(Grano);
  for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
  for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
  CalculeFuerzas(Grano);
  for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
  for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
  CalculeFuerzas(Grano);
  for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
  for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
  CalculeFuerzas(Grano);
  for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
  for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
}

//----------------- Funciones de Animacion ----------

class Animacion{
  private:
    double xmin,xmax,ymin,ymax;
    double Lx,Ly;
  public:
    void Init(double xmin0,double xmax0, double ymin0,double ymax0,
              double Lx0,double Ly0){
      xmin=xmin0 ;xmax=xmax0;
      ymin=ymin0; ymax=ymax0;
      Lx = Lx0; Ly=Ly0;
    }
    void InicieAnimacion(ofstream &file);
    void InicieCuadro(ofstream &file);
    void TermineCuadro(ofstream &file);

};

void Animacion::InicieAnimacion(ofstream &file){
  file<<"set terminal gif animate"<<endl;
  file<<"set output '2DGas.gif'"<<endl;
  file<<"unset key"<<endl;
  file<<"set xrange["<<xmin<<":"<<xmax<<"]"<<endl;
  file<<"set yrange["<<ymin<<":"<<ymax<<"]"<<endl;
  file<<"set size ratio -1"<<endl;
  file<<"set parametric"<<endl;
  file<<"set trange [0:7]"<<endl;
  file<<"set isosamples 12"<<endl;
}
void Animacion::InicieCuadro(ofstream &file){
    file<<"plot 0,0 ";
    file<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    file<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    file<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    file<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void Animacion::TermineCuadro(ofstream &file){
    file<<endl;
}
