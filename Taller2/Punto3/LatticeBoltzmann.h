#include <cmath>
#include <iostream>
#include <fstream>
#include "vector.h"

#ifndef __LATTICEBOLTZMANN_H_
#define __LATTICEBOLTZMANN_H_

using namespace std;

// Dimensiones del sistema
const int Lx=512;
const int Ly=64;

// Parámetros del LB
const int Q=9; //Direcciones en la celda

const double tau=1.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

const double nu=(tau-0.5)/3; // Viscosidad cinemática

//------------------ Clase LatticeBoltzmann ------------------
class LatticeBoltzmann{
private:
  double w[Q]; // Pesos
  int V[2][Q]; // Vectores de velocidad
  double *f, *fnew; // Functiones de distribución
  
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix, int iy, int i){return (ix*Ly*Q)+(iy*Q)+i;};
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double feq(double rho0, double Ux0, double Uy0, int i);
  void Start(double rho0, double Ux0, double Uy0);
  void Collision(void);
  void ImposeFields(double Ufan, int ixc, int iyc, int R);
  void Advection(void);
  void Print(const char * NameFile, double Ufan);
  double dUa_dXb(int ix, int iy, int a, int b);
  double sigmaxx(int ix, int iy);
  double sigmayy(int ix, int iy);
  double sigmaxy(int ix, int iy);
  double interp(double x, double y, int ix, int iy, int a);
  vector3D calcular_dF(double x, double y, vector3D dA);
  vector3D fuerzas(int N, int ixc, int iyc, int R);
  double print_deriv(const char * NameFile);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  // Asignar los pesos
  w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
  // Asignar vectores de velocidad
  V[0][0]=0; V[0][1]=1; V[0][2]=0; V[0][3]=-1; V[0][4]=0;
  V[1][0]=0; V[1][1]=0; V[1][2]=1; V[1][3]=0;  V[1][4]=-1;

  V[0][5]=1; V[0][6]=-1; V[0][7]=-1;  V[0][8]=1;
  V[1][5]=1; V[1][6]=1;  V[1][7]=-1;  V[1][8]=-1;
  // Crear arreglos dinámicos
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize]; fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f; delete[] fnew;
}
double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=V[0][i]*fnew[n0]; else sum+=V[0][i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=V[1][i]*fnew[n0]; else sum+=V[1][i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) // Para cada celda
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){
	n0=n(ix,iy,i);
	f[n0]=feq(rho0,Ux0,Uy0,i);
      }
}
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++) // Para cada celda
    for(iy=0;iy<Ly;iy++){
      // Calcule campos macroscópicos en cada celda
      rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){ // Para cada vector de velocidad
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0] + Utau*feq(rho0,Ux0,Uy0,i);
      }
    }
}
void LatticeBoltzmann::ImposeFields(double Ufan, int ixc, int iyc, int R){
  int i,ix,iy,n0; double rho0, R2=R*R;
  // Recorre todas las celdas, revisando si son abanico u obstáculo
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,false);
      // Abanico
      if(ix==0)
	for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,Ufan,0,i);}
      // Obstáculo
      else if((ix-ixc)*(ix-ixc)+(iy-iyc)*(iy-iyc)<=R2)
	for(i=0;i<Q;i++){n0=n(ix,iy,i); fnew[n0]=feq(rho0,0,0,i);}
    }
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ // En cada dirección
	ixnext=(ix+V[0][i]+Lx)%Lx; iynext=(iy+V[1][i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; // Fronteras períodicas
      }
}
void LatticeBoltzmann::Print(const char * NameFile, double Ufan){
  ofstream MyFile(NameFile); double rho0,Ux0,Uy0; int ix,iy;
  for(ix=0;ix<Lx;ix+=2){
    for(iy=0;iy<Ly;iy+=2){
      rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix,iy,true)/rho0;
      MyFile << ix << " " << iy << " " << Ux0/Ufan*4 << " " << Uy0/Ufan*4 << endl;
    }
    MyFile << endl;
    }
  MyFile.close();
}
double LatticeBoltzmann::dUa_dXb(int ix, int iy, int a, int b){
  /*
    Derivada de la componente * a * de U, con respecto
    a la coordenada *X_b*.

    Parámetros
    ----------
    ix : Índice X de la celda
    iy : Índice Y de la celda
    a : Componente de la velocidad. 0 -> X
                                    1 -> Y
    b : Coordenada con respecto a la cual se deriva. 0 -> X
                                                     1 -> Y
  */
  int ii, ix_next, iy_next;
  double rho0, sum, Ua;
  for(ii=0,sum=0; ii<Q; ii++){
    ix_next=ix+V[0][ii]; iy_next=iy+V[1][ii];
    rho0=rho(ix_next,iy_next,true);
    if(a==0) Ua=Jx(ix_next,iy_next,true)/rho0;
    else Ua=Jy(ix_next,iy_next,true)/rho0;
    sum+=w[ii]*V[b][ii]*Ua;
  }
  return 3*sum;
}
double LatticeBoltzmann::sigmaxx(int ix, int iy){
  double p, dUx_dx, rho0, eta;
  rho0=rho(ix,iy,true); p=rho0/3; eta=rho0*nu;
  dUx_dx=dUa_dXb(ix,iy,0,0);
  return -p+2*eta*dUx_dx;
}
double LatticeBoltzmann::sigmayy(int ix, int iy){
  double p, dUy_dy, rho0, eta;
  rho0=rho(ix,iy,true); p=rho0/3; eta=rho0*nu;
  dUy_dy=dUa_dXb(ix,iy,1,1);
  return -p+2*eta*dUy_dy;
}
double LatticeBoltzmann::sigmaxy(int ix, int iy){
  double dUx_dy, dUy_dx, rho0, eta;
  rho0=rho(ix,iy,true); eta=rho0*nu;
  dUx_dy=dUa_dXb(ix,iy,0,1); dUy_dx=dUa_dXb(ix,iy,1,0);
  return eta*(dUx_dy+dUy_dx);
}
double LatticeBoltzmann::interp(double x, double y, int ix, int iy, int a){
  /*
    Toma los valores del campo en las esquinas de una celda, e interpola
    al punto interior (x,y).
    La variable * a * señala cuál componente del tensor de esfuerzos se usa:
    a==0 -> sigmaxx
    a==1 -> sigmaxy
    a==2 -> sigmayy
  */
  double phi_00, phi_10, phi_01, phi_11, u, v;
  if(a==0){
    phi_00=sigmaxx(ix,iy); phi_01=sigmaxx(ix,iy+1);
    phi_10=sigmaxx(ix+1,iy); phi_11=sigmaxx(ix+1,iy+1);
  } else if(a==1){
    phi_00=sigmaxy(ix,iy); phi_01=sigmaxy(ix,iy+1);
    phi_10=sigmaxy(ix+1,iy); phi_11=sigmaxy(ix+1,iy+1);
  } else {
    phi_00=sigmayy(ix,iy); phi_01=sigmayy(ix,iy+1);
    phi_10=sigmayy(ix+1,iy); phi_11=sigmayy(ix+1,iy+1); 
  }
  u=(x-ix); v=(y-iy);
  return phi_00*(1-u)*(1-v)+phi_10*u*(1-v)+phi_01*(1-u)*v+phi_11*u*v;
}
vector3D LatticeBoltzmann::calcular_dF(double x, double y, vector3D dA){
  // Índices de la esquina inferior izquierda de la celda
  int ix=(int) x; int iy=(int) y;
  // Componentes del tensor de esfuerzo
  double sigma_00, sigma_01, sigma_11;
  sigma_00=interp(x,y,ix,iy,0);
  sigma_01=interp(x,y,ix,iy,1);
  sigma_11=interp(x,y,ix,iy,2);
  // Vectores columna del tensor de esfuerzo
  vector3D sigma_row_0; sigma_row_0.load(sigma_00,sigma_01,0);
  vector3D sigma_row_1; sigma_row_1.load(sigma_01,sigma_11,0);
  // Componentes del vector fuerza por unidad de longitud
  double Fx, Fy;
  Fx=sigma_row_0*dA; Fy=sigma_row_1*dA;
  // Vector fuerza por unidad de longitud
  vector3D dF; dF.load(Fx,Fy,0);
  return dF;
}
vector3D LatticeBoltzmann::fuerzas(int N, int ixc, int iyc, int R){
  /*
    Inscribimos un polígono regular de N lados en el obstáculo circular.
    
  */
  // Angulo que substiende cada segmento
  double theta=2*M_PI/N;
  // Longitud del lado
  double l=2*R*sin(theta/2);
  // Apotema
  double a=R*cos(theta/2);
  // Vector de área y de fuerza por unidad de longitud
  vector3D dA, dF;
  // Vector de fuerza resultante
  vector3D F; F.load(0, 0, 0);
  
  for(int i=0; i<N; i++){
    double x=ixc+a*cos(i*theta), y=iyc+a*sin(i*theta);
    dA.load(l*cos(i*theta), l*sin(i*theta), 0);
    dF=calcular_dF(x, y, dA);
    F+=dF;
  }

  return F;
}
double LatticeBoltzmann::print_deriv(const char * NameFile){
  ofstream MyFile(NameFile); int ix,iy;
  for(ix=140,iy=32; ix<Lx; ix++){
    MyFile << ix << " " << dUa_dXb(ix,iy,0,0) << endl;
  }
  MyFile.close();
}

#endif
