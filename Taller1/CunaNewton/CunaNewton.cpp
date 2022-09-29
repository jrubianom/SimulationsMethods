#include <iostream>
#include <cmath>
#include "EvolucionCuna.h"

//-------------Constantes-----------------
//Constantes fisica
const double g = 980;
const double alpha = 1.5;
//Numero esferas y tiempos de simiulacion
const int N = 3;
const int TUdt = 1e6; //Razon entre T y dt
const double NT = 0.28;//Numero de periodos que correrá la simulacion
const int TUtcuadro = 10; //Razon entre Periodo y tCuadro
const bool gif = false; //Set true solo si corre para un solo K
//Propiedades Esferas
const double m = 100,l=12,R=1.5,theta0=-15*M_PI/180.0;
const int ks = 6;
const double K[ks] = {0.1e10,0.2e10,0.5e10,1e10,2.5e10,10e10};

//------------------------------------------
using namespace std;

int main(){

  //Archivo que almacena datos de t, tau
  string file = "Exp"+to_string(0)+".txt";
  //----Parametros--------
  Parametros P;
  P.Inicie(m,l,R,theta0,K[0],N,TUdt,TUtcuadro);
  P.g = g; P.alpha = alpha; P.NT = NT;
  P.gif = gif;
  //----Simulacion------------
  //Pendulo al que se le medirá el torque
  int penduloID = 1;

  for(int i=0; i < ks; i++){
    P.K = K[i];
    Evolucion(penduloID,P,file);
    //Crear nuevo file para la actual K
    file = "Exp"+to_string(i+1)+".txt";
  }

  return 0;
}
