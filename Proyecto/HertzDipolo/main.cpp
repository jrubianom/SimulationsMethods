#include<iostream>
#include<cmath>
#include <fstream>
#include "LBEM.h"

using namespace std;

//------------------------CONSTANTS-------------------------------
const int Lx = 100;   //
const int Ly = 100;   //
const int Lz = 100; //
const int Qr = 2, Qp = 3, Qi = 4, Qj = 2;
//-------------------
const double Tau = 0.5;
const double UTau = 1/Tau;
const double UmUTau=1-1/Tau;
//-------------------
const double Epsilon0=1, Mu0=2;
const double Sigma0=0.0;
const double C=1.0/sqrt(2.0);

const double E00=0.001,B00=E00/C;
const double J0=0.0001;

const double alpha=0.75;
const double T=25;

int main(){

  Parameter Params(Lx,Ly,Lz,Tau,Epsilon0,Mu0,Sigma0,E00,B00,J0,alpha,T);
  LatticeBoltzmann Dipole(Params);
  int t, tmax=70;
  
  Dipole.Start();
  
  for(t=0;t<tmax;t++){
    cout<< "Iteracion:\t"<< t+1 << endl;
    Dipole.Collision(t);
    Dipole.Advection();
  }
  
  Dipole.Print();
  
  return 0;
}