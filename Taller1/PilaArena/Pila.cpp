// Simular el movimiento de N granos con gravedad, fricción y choques inelásticos
#include <iostream>
#include <cmath>
#include "/home/live/repos/Files/vector.h"
#include "Random64.h"
#include "fstream"
#include <string>

using namespace std;

//---- declarar constantes ---
const double g=9.8, K=1.0e4, Gamma=150, Kcundall=500, MU=0.4;
const double Lx=160, Ly=100, L0 = 60;
const int N=350;
const int Nsadd = 0,Nss = 80,Ns = Nss+Nsadd,  Nt = Ns + N + 2;
const double Rs = Lx/(2*Nss),Rr = 2;

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
    vector3D r,V,F; double m,R; double theta,omega,tau; double I;
  public:
    void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
                double theta0,double omega0);
    void BorreFuerza(){F.load(0,0,0); tau=0;};
    void AdicioneFuerza(vector3D F0){F+=F0;};
    void AdicioneTorque(double tau0){tau+=tau0;};
    void Mueva_r(double dt, double Coeficiente);
    void Mueva_V(double dt, double Coeficiente);
    void Dibujese(void);
    double Getx(void){return r.x();}; //inline
    double Gety(void){return r.y();}; //inline
    double Gettheta(void){return theta;}; //inline
    friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
                    double theta0,double omega0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
}
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);  omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";
}

//--- clase Colisionador ----
class Colisionador{
  private:
    double dcontacto[Nt][Nt],hold[Nt][Nt];
    int Nlevel = N;
  public:
    void Inicie();
    void CalculeFuerzas(Cuerpo * Grano,double dt);
    void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,
                            double &dcontacto,double &hold,double dt);
    void SetLevel(int nlevel){Nlevel = nlevel;}
    int Getlevel(){return Nlevel;}
};

void Colisionador::Inicie(){
  for(int i=0;i<N;i++){
    for(int j=0;j<Nt;j++){
      dcontacto[i][j] = hold[i][j] = 0;
    }
  }
}

void Colisionador::CalculeFuerzas(Cuerpo * Grano, double dt){
  int i,j; vector3D Fg;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<Nlevel;i++)
    Grano[i].BorreFuerza();
  for(i=N;i<Nt;i++)
    Grano[i].BorreFuerza();
  //--- Sumar el peso ---
  for(i=0;i<Nlevel;i++){
    Fg.load(0,-Grano[i].m*g,0);
    Grano[i].AdicioneFuerza(Fg);
  }
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<Nlevel;i++){
    for(j=i+1;j<Nlevel;j++)
      CalculeFuerzaEntre(Grano[i], Grano[j],dcontacto[i][j],hold[i][j],dt);
    for(j=N;j<Nt;j++)
      CalculeFuerzaEntre(Grano[i], Grano[j],dcontacto[i][j],hold[i][j],dt);
  }
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2,
                                      double &dcontacto, double &hold,double dt){
  vector3D r21=Grano2.r-Grano1.r;
  double R1 = Grano1.R,R2 = Grano2.R;
  double d=r21.norm(), s=R1+R2-d;
  if(s>0){
    double m1 = Grano1.m, m2 = Grano2.m, m12 = m1*m2/(m1+m2);
    vector3D n=r21*(1.0/d);
    vector3D V21= Grano2.V - Grano1.V;
    vector3D Rw; Rw.load(0,0,R1*Grano1.omega + R2*Grano2.omega);
    vector3D Vc = V21-(Rw^n),t;
    t.load(n.y(),-n.x(),0);
    double Vn = Vc*n, Vt = Vc*t;
    //Fuerzas centrales
    double Fn = ((K*pow(s,1.5))-m12*sqrt(s)*Gamma*Vn);
    vector3D F2=n*Fn;
    Grano2.AdicioneFuerza(F2);   Grano1.AdicioneFuerza(F2*(-1));
    dcontacto+=Vt*dt;
    double Ft=-Kcundall*dcontacto;
    double Ftmax = MU*fabs(Fn);
    if(fabs(Ft)> Ftmax) Ft = Ft/fabs(Ft)*Ftmax;

    Grano2.AdicioneFuerza(t*Ft); Grano1.AdicioneFuerza(t*(Ft*(-1)));
    Grano2.AdicioneTorque(R2*Ft); Grano1.AdicioneTorque(R1*Ft);
  }

  if(hold>0 && s<0) dcontacto = 0;
  hold = s;
}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl;
  cout<<"set output 'SandPile.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange["<<-Nsadd/2.0 <<":"<<Lx+Nsadd/2.0<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
  cout<<"plot 0,0 ";
  cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
  //cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
  cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
  cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
  cout<<endl;
}

//-----------Funciones Para la evolucion de un paso --------------------
void EvolucionUnGrano(double tmax,double dt,double tcuadro,
                      Cuerpo *Grano, Colisionador &Hertz);

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Grano[Nt];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1, R0=Rr, kT=0;
  double tmax=5*sqrt(Ly/g),tcuadro=tmax/10,dt=1e-2;
  double Omega, omega_max = 8.0;//sqrt(Ly*g)/R0*1e-9;
  
  InicieAnimacion(); //Dibujar

  //Inicializar las paredes
  double Rpared=100*Lx, Mpared=100*m0;
  //---------------(  x0,       y0,Vx0,Vy0,    m0,    R0, theta0,omega0)
  Grano[N+Ns].Inicie(Lx+Rpared,Ly/2,  0,  0,Mpared,Rpared,      0,     0); //Pared derecha
  Grano[N+Ns+1].Inicie(  -Rpared,Ly/2,  0,  0,Mpared,Rpared,      0,     0); //Pared izquierda
  //Grano[N+Ns+2].Inicie(Lx/2,Ly+Rpared,  0,  0,Mpared,Rpared,      0,     0); //Pared de arriba
  for(int l=N, ix=0; l < N+Ns; l++,ix++)
    Grano[l].Inicie((2*ix+1-Nsadd)*Rs,  -Rs,  0,  0,Mpared,Rs,      0,     0); //Pared de abajo
  //Inicializar las moléculas
  for(int i=0;i<N;i++){
      Omega=omega_max*ran64.r();
      //--------------------(   x0,   y0,          Vx0,          Vy0, m0,R0,theta0,omega0)
      Grano[i].Inicie(Lx/2,L0-2*R0,0,0, m0,R0,0,Omega);//OJO
  }

  ofstream file("Nivel.txt");
  for(int Nlevel = 1;Nlevel<= N; Nlevel++){
    Hertz.SetLevel(Nlevel);
    EvolucionUnGrano(tmax,dt,tcuadro,Grano,Hertz);
    file << Nlevel << endl;
  }
  file.close();

  return 0;
}

void EvolucionUnGrano(double tmax,double dt,double tcuadro,
                      Cuerpo *Grano, Colisionador &Hertz){
  double tdibujo,t;
  int i,Nlevel = Hertz.Getlevel();
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){
    //Dibujar
    if(tdibujo>tcuadro){
      InicieCuadro();
      for(i=0;i<Nlevel;i++) Grano[i].Dibujese();
      for(i=N;i<N+Ns;i++) Grano[i].Dibujese();
      TermineCuadro();
      tdibujo=0;

      if(Nlevel == N){
        cout<<"set term pdf"<<endl;
        cout<<"set output 'Final.pdf'"<<endl;
        cout<<"unset key"<<endl;
        cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
        cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
        cout<<"set size ratio -1"<<endl;
        cout<<"set parametric"<<endl;
        cout<<"set trange [0:7]"<<endl;
        cout<<"set isosamples 12"<<endl;

        InicieCuadro();
        for(i=0;i<N+Ns;i++) Grano[i].Dibujese();
        TermineCuadro();
      }
    }

    //--- Muevase por PEFRL ---
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano,dt);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<Nlevel;i++)Grano[i].Mueva_r(dt,epsilon);
  }
}
