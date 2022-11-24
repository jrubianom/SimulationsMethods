#include <cmath>
#include <iostream>
#include <fstream>
#include "LatticeBoltzmann.h"

using namespace std;

/*
  El siguiente programa calcula la fuerza de Magnus para varias velocidades 
  del fluido (Ufan0). En cada iteración computa también el número de Reynolds e imprime
  los resultados a un archivo de datos ("plot_data.dat")
*/

// Posición y tamaño del cilindro obstáculo
const int ixc=Lx/4, iyc=Ly/2, R=Ly/8;
const double R2=R*R;

//------------------   Funciones globales   ------------------
double ReynoldsNumber(double U, double L, double k_viscosity){
  return U*L/k_viscosity;
}
double fuerzaMagnus(double rho0, double U0, double omega0){
  return 0.5*rho0*2*R*R*omega0*U0;
}
//---------------- Función main ------------------ 
int main(void){
  
  int t, tmax=10000;
  double rho0=1.0, omega0=2*M_PI/1000, Ufan0;

  ofstream data_file("plot_data.dat");
    
  
  for(Ufan0=0.10; Ufan0<0.30; Ufan0+=0.01){
    LatticeBoltzmann Aire;
    Aire.Start(rho0, Ufan0, 0);
    // Correr
    for(t=0;t<tmax;t++){
      Aire.Collision();
      Aire.ImposeFields(Ufan0, omega0, ixc, iyc, R);
      Aire.Advection();
    }
    // Calcular fuerza de arrastre
    vector3D F; F=Aire.fuerzas(24, ixc, iyc, R);
    double Rn=ReynoldsNumber(U, 2*R, nu);
    double Flift=F.y();
    double Fdrag=F.x();
    // Imprimir
    data_file << Ufan0 << " " << Rn << " " << Flift << endl;
  }
  data_file.close();
  
  return 0;
}
