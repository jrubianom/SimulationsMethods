#include <iostream>
#include <cmath>
#include "EvolucionCuna.h"
#include <algorithm>
#include <iterator>

//-------------COnstantes-----------------
//Constantes fisica
const double g = 980;
const double alpha = 1.5;
//Numero esferas y tiempos de simiulacion
const int N = 3;
const int TUdt = 1e5; //Razon entre T y dt
const double NT = 0.28;//1.0/4.0; //Numero de periodos que correrá la simulacion
const int TUtcuadro = 10; //Razon entre Periodo y tCuadro
const bool gif = false;
//Propiedades Esferas
const double m = 100,l=12,R=1.5,theta0=-15*M_PI/180.0;
const int ks = 6;
const double K[ks] = {0.1e10,0.2e10,0.5e10,1e10,2.5e10,10e10};

//------------------------------------------
using namespace std;

void printMax(double *Torques,Parametros P,int N_iter){
  double *maxtau = max_element(Torques, Torques+N_iter);
  double *mintau = min_element(Torques, Torques+N_iter);
  int t_maxtau = distance(Torques,maxtau);
  int t_mintau = distance(Torques,mintau);
  double tmax = P.dt*(t_mintau-t_maxtau);
  cout << P.K <<" " << *maxtau << " "<< tmax<< endl;
}

int main(){

  string file = "Exp"+to_string(0)+".txt";
  //----Parametros--------
  Parametros P;
  P.Inicie(m,l,R,theta0,K[0],N,TUdt,TUtcuadro);
  P.g = g; P.alpha = alpha; P.NT = NT;
  P.gif = gif;
  //----Simulacion------------
  int N_iter = (int) (NT*TUdt);
  double Torques[N_iter];
  int penduloID = 1;

  for(int i=0; i < ks; i++){
    Evolucion(Torques,penduloID,P,file);
    P.K = K[i];
    //-------Imprimir taus maximos y tmax---------
    printMax(Torques,P,N_iter);
    //actualizar nombre de file
    file = "Exp"+to_string(i+1)+".txt";
  }

  return 0;
}
