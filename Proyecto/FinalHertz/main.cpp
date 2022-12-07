#include <algorithm>
#include<iostream>
#include<cmath>
#include <fstream>
#include <vector>
#include <functional>
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

const double alpha=0.5;
const double T=25;
const double lamd = C*T;
const double Z0=sqrt(Mu0/Epsilon0);

void PrintPower(vector<double> &V,int N,double R,string name,Theo &formulas,bool flag);
void compare(vector<double> &V1,vector<double> &V2,int N);

void InicieCuadro(ofstream &file,string frame){
  file<<"plot '" << frame <<"'"<< " u 1:3 w l title 'teorico','' u 1:2 title 'LB'"<<endl;
}

void TermineCuadro(ofstream &file){
    file<<endl;
}

void InicieAnimacionB(ofstream &file){
  file<<"set terminal gif animate"<<endl;
  file<<"set output 'FiledBy.gif'"<<endl;
  file<<"unset key"<<endl;
  file << "set g" <<endl;
  file<<"set xrange[55:"<<Lx<<"]"<<endl;
  file << "set yrange[-0.2:0.2]" <<endl;
}

void InicieAnimacionE(ofstream &file){
  file<<"set terminal gif animate"<<endl;
  file<<"set output 'FiledEz.gif'"<<endl;
  file<<"unset key"<<endl;
  file<<"set xrange[55:"<<Lx<<"]"<<endl;
  file << "set yrange[-0.1:0.1]" <<endl;
}

void InicieAnimacionS(ofstream &file){
  file<<"set terminal gif animate"<<endl;
  file<<"set output 'FiledS.gif'"<<endl;
  file<<"unset key"<<endl;
  file<<"set xrange[55:"<<Lx<<"]"<<endl;
  file << "set yrange[0:0.2]" <<endl;
}


int main(){

  Parameter Params(Lx,Ly,Lz,Tau,Epsilon0,Mu0,Sigma0,E00,B00,J0,alpha,T);
  LatticeBoltzmann Dipole(Params);
  int t, tmax=3*T;
  int N=100;
  vector<double> S(N,0),S_current(N,0),Splus(N,0);
  vector<double> S2(N,0),S2_current(N,0),S2plus(N,0);
  double R = lamd;
  Dipole.Start();

  string current_frameB = "Datos/frameB.dat";
  ofstream AnimFileB("animationBy.gp");
  string current_frameE = "Datos/frameE.dat";
  ofstream AnimFileE("animationEz.gp");

  string current_frameS = "Datos/frameS.dat";
  ofstream AnimFileS("animationS.gp");


  InicieAnimacionB(AnimFileB); InicieAnimacionE(AnimFileE); InicieAnimacionS(AnimFileS);


  for(t=0;t<tmax;t++){
    Dipole.UpdateTime(t);
    cout<< "Iteracion:\t"<< t+1 <<endl;
    Dipole.Collision();
    Dipole.Advection();
    //Dipole.ImposeFields();

    if(R <= C*t && C*t <= R+lamd){
      //Dipole.PowerPlanePhi(R,N,S_current);
      //transform(S.begin(),S.end(),S_current.begin(),Splus.begin(),plus<double>());
      //S.assign(Splus.begin(),Splus.end());
      Dipole.PowerPlane(R,N,0,S_current);
      //transform(S.begin(),S.end(),S_current.begin(),Splus.begin(),plus<double>());
      //S.assign(Splus.begin(),Splus.end());
      compare(S,S_current,N);

      //Dipole.PowerPlane(R,N,M_PI/2,S2_current);
      Dipole.PowerPlaneTheta(R,N,S2_current);
      //transform(S2.begin(),S2.end(),S2_current.begin(),S2plus.begin(),plus<double>());
      //S2.assign(S2plus.begin(),S2plus.end());
      compare(S2,S2_current,N);


      Dipole.Print(current_frameB,current_frameE,current_frameS);
      InicieCuadro(AnimFileB,current_frameB);
      current_frameB = "Datos/frameB"+to_string(t)+".dat";

      InicieCuadro(AnimFileE,current_frameE);
      current_frameE = "Datos/frameE"+to_string(t)+".dat";

      InicieCuadro(AnimFileS,current_frameS);
      current_frameS = "Datos/frameS"+to_string(t)+".dat";

    }

  }

  //Dipole.Print();
  //Print Power data
  Theo formulas; formulas.Init(Params);
  PrintPower(S,N,R,"EPlane.txt",formulas,1);
  PrintPower(S2,N,R,"BPlane.txt",formulas,0);

  return 0;
}

void PrintPower(vector<double> &V,int N,double R,string name,Theo &formulas,bool flag){
  ofstream file("Datos/"+name);
  ofstream file2("Datos/Teo"+name);
  double dphi,phi;
  double Pteo;
  dphi = 2*M_PI/(N-1);
  for(int i=0; i < N; i++){
    phi = i*dphi;
    if(flag){
      Pteo = formulas.Power(phi,0,R);
      file2 << M_PI/2-phi  << " " << Pteo/(Z0*J0*J0) << endl;
    }
    else{
      Pteo = formulas.Power(M_PI/2,phi,R);
      file2 << phi << " " << Pteo/(Z0*J0*J0) << endl;
    }
    file << M_PI/2-phi << " " << V[i]/(Z0*J0*J0)*0.5 << endl;
  }
  file.close();
  file2.close();
}


void compare(vector<double> &V1,vector<double> &V2,int N){
  for(int i=0; i < N; i++)
    if(V1[i] < V2[i])
      V1[i] = V2[i];
}
