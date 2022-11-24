#include <cmath>
#include <iostream>
#include <fstream>
#include "LatticeBoltzmann.h"

using namespace std;

/*
  El siguiente programa calcula la fuerza de arrastre y el coeficiente
  de arrastre para varias   velocidades del fluido (Ufan0).
  Además, en cada iteración computa también el número de Reynolds e imprime
  los resultados a un archivo de datos ("dragCoeff_vs_reynosldNumber.dat")
  Primera columna: Número de Reynolds
  Segunda columna: Coeficiente de arrastre
  Se ejecuta el archivo .py para producir la gr
*/

// Posición y tamaño del cilindro obstáculo
const int ixc=Lx/4, iyc=Ly/2, R=Ly/8;
const double R2=R*R;

//------------------   Funciones globales   ------------------
double ReynoldsNumber(double U, double L, double k_viscosity){
  return U*L/k_viscosity;
}

//---------------- Función main ------------------ 
int main(void){
  
  int t, tmax=10000;
  double rho0=1.0, Ufan0;

  ofstream data_file("dragCoeff_vs_reynosldNumber.dat");
  
  for(Ufan0=0.1; Ufan0<0.3; Ufan0+=0.01){
    LatticeBoltzmann Aire;
    Aire.Start(rho0, Ufan0, 0);
    // Correr
    for(t=0;t<tmax;t++){
      Aire.Collision();
      Aire.ImposeFields(Ufan0, ixc, iyc, R);
      Aire.Advection();
    }
    double R_n=ReynoldsNumber(Ufan0, 2*R, nu);
    // Calcular fuerza de arrastre
    vector3D F; F=Aire.fuerzas(24, ixc, iyc, R);
    double Fdrag=F.x();
    // Coeficiente de arrastre
    double C=2*Fdrag/(rho0*2*R*Ufan0);
    double C2=2*Fdrag/(rho0*2*R*Ufan0*Ufan0);
    // Imprimir
    data_file << R_n << " " << C2 << endl;
  }
  data_file.close();
  
  return 0;
}
