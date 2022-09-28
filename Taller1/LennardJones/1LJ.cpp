#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
#include "fstream"

using namespace std;

//---- declarar constantes ---
const double g=0, eps=1.0, r0 = 10, R = 2.5;
const int N=1;
const bool gif = true;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
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
    double Energia(void);
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
  file<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

double Cuerpo::Energia(void){
  double E = 0;
  double rr = pow(r0/r.norm(),6);
  E += 0.5*m*V.norm2();
  E += eps*rr*(rr - 2);
  return E*1e3;
}

//--- clase Colisionador ----
class Colisionador{
  private:
  public:
    void CalculeFuerzas(Cuerpo &Particula);
    void CalculeFuerzaCentral(Cuerpo & Particula);
};

void Colisionador::CalculeFuerzas(Cuerpo &Particula){
  vector3D Fg;
  //--- Borrar todas las fuerzas ---
  Particula.BorreFuerza();
  //--- Sumar el peso ---
  Fg.load(0,-Particula.m*g,0);
  Particula.AdicioneFuerza(Fg);
//--- Calcular Fuerzas entre pares de granos ---
  CalculeFuerzaCentral(Particula);
}

void Colisionador::CalculeFuerzaCentral(Cuerpo & Particula1){
  vector3D r=Particula1.r;
  double d=r.norm();
  vector3D n=r*(1.0/d);
  vector3D F=n*(12*eps/d)*(pow(r0/d,12)-pow(r0/d,6));
  Particula1.AdicioneFuerza(F);

}

//----------------- Funciones de Animacion ----------
void InicieCuadro(ofstream &file){
  file<<"plot 0,0 ";
}

void TermineCuadro(ofstream &file){
  file<<endl;
}

void InicieAnimacion(ofstream &file){

  file<<"set terminal gif animate"<<endl;
  file<<"set output 'UnaParticula.gif'"<<endl;
  file <<"set g" << endl;
  file<<"unset key"<<endl;
  file<<"set xrange[0:"<< 2*r0 << "]"<<endl;
  file<<"set yrange["<<-2*R<<":"<< 2*R<<"]"<<endl;
  file<<"set size ratio -1"<<endl;
  file<<"set parametric"<<endl;
  file<<"set trange [0:7]"<<endl;
  file<<"set isosamples 12"<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Particula;
  Colisionador LenJones;
  Crandom ran64(1);
  double m0=1, R0=R, kT=0.5, V0=sqrt(2*kT/m0);
  double x0 = 10, y0 = 0,vx0 = V0, vy0 = 0;
  int i;
  double t,tdibujo,tmax=100,tcuadro=tmax/1000,dt=1e-1;

  ofstream filegif("Animacion1.txt");
  ofstream filedata("datos1.txt");

  if(gif) InicieAnimacion(filegif); //Dibujar

  //Inicializar las moléculas
//--------------------(x0,y0, Vx0, Vy0, m0,R0,theta0,omega0)
  Particula.Inicie(x0,y0,vx0,vy0, m0,R0);

  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(gif && tdibujo>tcuadro){
      
      InicieCuadro(filegif);
      for(i=0;i<N;i++) Particula.Dibujese(filegif);
      TermineCuadro(filegif);
      
      tdibujo=0;
    }


    filedata << t << "\t"<< Particula.Getx() <<"\t" << Particula.Energia()<<endl;

    //--- Muevase por PEFRL ---
    Particula.Mueva_r(dt,epsilon);
    LenJones.CalculeFuerzas(Particula);
    Particula.Mueva_V(dt,lambda2);
    Particula.Mueva_r(dt,chi);
    LenJones.CalculeFuerzas(Particula);
    Particula.Mueva_V(dt,lambda);
    Particula.Mueva_r(dt,chiepsilon);
    LenJones.CalculeFuerzas(Particula);
    Particula.Mueva_V(dt,lambda);
    Particula.Mueva_r(dt,chi);
    LenJones.CalculeFuerzas(Particula);
    Particula.Mueva_V(dt,lambda2);
    Particula.Mueva_r(dt,epsilon);

  }
  filegif.close();
  filedata.close();

  return 0;
}
