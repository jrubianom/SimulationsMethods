#include <cstdint>
#include <fstream>
#include <iostream>
#include <cmath>
#include "vector.h"
#include "fstream"
#include <string>

using namespace std;

//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
  public:
    //vector3D r,V,F;  double m,R;
    double theta,omega,tau;
    double m,l,R,x0,I;
  public:
    void Inicie(double theta0, double omega0, double m0,
                double R0,double l0, double x0);
    void CalculeMomentoI(){I = l*l*m;};
    void BorreTorque(void){tau = 0;};
    void SumeTorque(double tau0){tau+=tau0;};
    void Mueva_theta(double dt,double coeficiente);
    void Mueva_omega(double dt,double coeficiente);
    void Dibujese(ofstream &file);
    double Getx(void){return x0+l*sin(theta);}; //Inline
    double Gety(void){return -l*cos(theta);}; //Inline
    double Gettau(void){return tau;}; //Inline
    friend class Colisionador;
};
void Cuerpo::Inicie(double theta0, double omega0, double m0,
                    double R0,double l0, double x00){
  theta = theta0;
  omega = omega0;
  m=m0; R=R0; l = l0; x0 = x00;
  CalculeMomentoI();
}
void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta += omega*(dt*coeficiente);
}
void Cuerpo::Mueva_omega(double dt,double coeficiente){
  omega+=tau*(dt*coeficiente/I);
}

void Cuerpo::Dibujese(ofstream &file){
  file<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
  file<<" , "<<x0<<"+"<<l<<"*t/7*sin("<<theta<<"),"<<l<<"*(-1)*t/7*cos("<<theta<<")";
}

//---------- Clase Colisionador --------------
class Colisionador{
  private:
    double g,K,alpha;
    int N;
  public:
    void Inicie(double g0,double K0,double alpha0,int N0){
      g=g0;K=K0;alpha=alpha0;
      N = N0;
    }
    void CalculeTorques(Cuerpo * Pendulo);
    void CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);
};
void Colisionador::CalculeTorques(Cuerpo * Pendulo){
  int i;
  double tau0;
  //Borrar fuerzas
  for(i=0;i<N;i++){
    Pendulo[i].BorreTorque();
    tau0 = -Pendulo[i].l*Pendulo[i].m*g*sin(Pendulo[i].theta);
    Pendulo[i].SumeTorque(tau0);
  }

  //Calcular torques entre todas las parejas de pendulos
  for(i=0;i<N-1;i++)
    CalculeTorqueEntre(Pendulo[i],Pendulo[i+1]);
}

void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  vector3D r1,r2,r,F,Tau1,Tau2,x;
  x.load(1,0,0);
  double s,rn;
  r1.load(Pendulo1.Getx(),Pendulo1.Gety(),0);
  r2.load(Pendulo2.Getx(),Pendulo2.Gety(),0);
  r = r2-r1; rn = r.norm();
  s = Pendulo1.R + Pendulo2.R - r.norm();
  if(s>0){
    F= r*(K*pow(s,alpha)/rn);
    r1-= (Pendulo1.x0*x);
    r2-= (Pendulo2.x0*x);
    Tau1 = F^r1; Tau2 = r2^F;
    Pendulo1.SumeTorque(Tau1.z()); Pendulo2.SumeTorque(Tau2.z());
  }
}

//------------------------------Clase Parametros -----------------
class Parametros{
  public:
    //Propiedades Esferas
    double m,l,R,theta0,K;
    double dt,T;
    //Constantes fisicas
    double g = 980;
    double alpha = 1.5;
    //Numero esferas y tiempos de simiulacion
    int N;
    int TUdt; //Razon entre T y dt
    double NT = 1; //Numero de periodos que correrá la simulacion
    int TUtcuadro; //Razon entre Periodo y tCuadro
    bool gif = true;
    void Inicie(double m0,double l0, double R0,double theta00,
                double K0, int N0,int TUdt0,int TUtcuadro0){
      m=m0;l=l0;R=R0;theta0=theta00;K=K0;N=N0;
      TUdt = TUdt0;
      TUtcuadro = TUtcuadro0;
    }
};

//----------- Funciones Animación -----------
void InicieCuadro(ofstream &file){
  file<<"plot 0,0 ";
}

void TermineCuadro(ofstream &file){
  file<<endl;
}

void InicieAnimacion(ofstream &file){

  file<<"set terminal gif animate"<<endl;
  file<<"set output 'Pendulos.gif'"<<endl;
  file<<"unset key"<<endl;
  file<<"set xrange[-10:30]"<<endl;
  file<<"set yrange[-20:0]"<<endl;
  file<<"set size ratio -1"<<endl;
  file<<"set parametric"<<endl;
  file<<"set trange [0:7]"<<endl;
  file<<"set isosamples 12"<<endl;
}

//--Funcion que ayuda a guardar los torques maximos, minimos y los respectivos timepos
void CantidadesMaximas(double &taumax,double &taumin,double
                       &timemax,double &timemin,
                       double tau, double t);

//------------------Evolucion del sistema--------------------------
void Evolucion(int penduloID,Parametros &P,string  fname){
  //Constantes
  const int N = P.N,TUdt = P.TUdt;
  const bool gif = P.gif;
  const double g=P.g, K=P.K,alpha=P.alpha;
  const double NT=P.NT;
  const int TUtcuadro = P.TUtcuadro;
  ///////////////////////////////
  Cuerpo Pendulo[N];
  Colisionador Newton; Newton.Inicie(g,K,alpha,N);
  double m0=P.m, l0=P.l, R=P.R, theta0 = P.theta0;
  double omega0=sqrt(g/l0), T=2*M_PI/omega0;
  double t, tmax = NT*T,dt=T/TUdt;
  P.dt = dt; P.T = T;
  double tdibujo,tcuadro = T/TUtcuadro;
  int i,j,jmax = (int) (NT*TUdt);
  //------Archivo para guardar Taus y tiempos
  ofstream file(fname);
  double tau,taumax,taumin,timemin,timemax;
  //---------------Inicializacion de los Pendulos------------------
  Pendulo[0].Inicie(theta0, 0, m0,  R,  l0, 0);
  for(i=1;i<N;i++){
    Pendulo[i].Inicie(0,0,m0,R,l0,2*i*R);
  }

  //------Guardar la animacion------
  string animaname= "AnimacionData.txt";
  ofstream AnimacionFile(animaname);
  if(gif){
    //Archivo para guardar datos de la animacion
    InicieAnimacion(AnimacionFile);
  }

//Evolucion del sistema
  for(t=0,tdibujo=0,j=0; j < jmax; t+=dt,tdibujo+=dt,j++){
//Dibujar
    if(tdibujo>tcuadro && gif){
      InicieCuadro(AnimacionFile);
      for(i=0;i<N;i++) Pendulo[i].Dibujese(AnimacionFile);
      TermineCuadro(AnimacionFile);
      tdibujo=0;
    }

    // Mover por PEFRL
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Coeficiente2);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);

    tau = Pendulo[penduloID].Gettau();
    if(t >= T/4.0) file << t <<" "<<tau<<" "<<endl;
    if(j == 0){
      taumax = tau;
      taumin = tau;
    }
    CantidadesMaximas(taumax,taumin,timemax,timemin,tau,t);
  }
  file.close();
  AnimacionFile.close();
  cout << K << " " << taumax << " " << (timemin-timemax)/2 << endl;
}

void CantidadesMaximas(double &taumax,double &taumin,double
                       &timemax,double &timemin,
                       double tau, double t){
  if(tau > taumax){
    taumax = tau;
    timemax = t;
  }
  if(tau <= taumin){
    taumin = tau;
    timemin = t;
  }
}
