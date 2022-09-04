#include <iostream>
#include <cmath>
#include "EvolucionCuna.h"
#include <algorithm>
#include <iterator>

using namespace std;

void FillTiempos(double *Tiempos,double dt,int N){
  for(int i = 0; i<N; i++){
    Tiempos[i] = dt*i;
  }
}

void print(double *Tiempos, double *Torques, int N){
  for(int i = 0; i < N; i++){
    cout << Tiempos[i] << " " << Torques[i] << endl;
  }
}

void printMax(double *Torques,Parametros P,int N_iter){
  double *maxtau = max_element(Torques, Torques+N_iter);
  double *mintau = min_element(Torques, Torques+N_iter);
  int t_maxtau = distance(Torques,maxtau);
  int t_mintau = distance(Torques,mintau);
  double tmax = P.dt*(t_mintau-t_maxtau);
  cout << P.K <<" " << *maxtau << " "<< tmax<< endl;
}

int main(){
  //Constantes
  //Propiedades Esferas
  double m = 100,l=12,R=1.5,theta0=-15*M_PI/180.0;
  int ks = 6;
  double K[ks] = {0.1e10,0.2e10,0.5e10,1e10,2.5e10,10e10};
  //Constantes fisica
  double g = 980;
  double alpha = 1.5;
  //Numero esferas y tiempos de simiulacion
  int N = 3;
  int TUdt = 1e5; //Razon entre T y dt
  double NT = 0.28;//1.0/4.0; //Numero de periodos que correrá la simulacion
  int TUtcuadro = 10; //Razon entre Periodo y tCuadro
  bool gif = false;
  string file = "Exp"+to_string(0)+".txt";
  //----Parametros--------
  Parametros P;
  P.Inicie(m,l,R,theta0,K[0],N,TUdt,TUtcuadro);
  P.g = g; P.alpha = alpha; P.NT = NT;
  P.gif = gif;
  //----Simulacion
  int N_iter = (int) (NT*TUdt);
  double Torques[N_iter],Tiempos[N_iter];
  int penduloID = 1;
  double *maxtau,*mintau;
  double tmax;
  int t_maxtau,t_mintau;
  for(int i=0; i < ks; i++){
    Evolucion(Torques,penduloID,P,file);
    P.K = K[i];
    if( i == 0 ) FillTiempos(Tiempos,P.dt,N_iter);
    //-------Imprimir taus maximos y tmax---------
    printMax(Torques,P,N_iter);
    //actualizar nombre de file
    file = "Exp"+to_string(i+1)+".txt";
  }
  //Guardar tmax,tmax y Plotear
  return 0;
}
