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
const double lamd = 16;
const double T=lamd/C;
const double Z0=sqrt(Mu0/Epsilon0);

const bool D3d = true;


////////////////////-----------/////////////
void PrintPower(vector<double> &V,int N,double R,string name,Theo &formulas,bool flag);

int main(){

  int ixs=Lx/2,iys=Ly/2,izs=Lz/2;
  Parameter Params(Lx,Ly,Lz,Tau,Epsilon0,Mu0,Sigma0,E00,B00,J0,alpha,T,
                   ixs,iys,izs);

  //Driven elements
  double lxs=0,lys=0,lzs=lamd*5.0/4.0,betaS=5;
  vector<Sources> JSources;
  JSources.push_back(Sources(ixs,iys,izs,lxs,lys,lzs,lamd,C,J0,betaS,alpha,D3d));

  //Parasitic elements
  double lxd=2,lyd=2,lzd=lamd*0.55; double lxr=2,lyr=2,lzr=lamd*0.56;
  double SigmaP = 0,beta = 0.5;
  double xd = lamd*0.125+ixs,spacing=0.2*lamd; double xr = ixs - lamd*0.35;
  vector<Parasitic> Paras;
  Paras.push_back(Parasitic(xd,iys,izs,lxd,lyd,0.55*lamd,SigmaP,beta)); //Director 1
  //Paras.push_back(Parasitic(xd+spacing,iys,izs,lxd,lyd,0.40*lamd,SigmaP,beta)); //Director 2
  //Paras.push_back(Parasitic(xd+2*spacing,iys,izs,lxd,lyd,0.35*lamd,SigmaP,beta)); //Director 3
  //Paras.push_back(Parasitic(xr,iys,izs,lxr,lyr,lzr,SigmaP,beta)); //Reflector


  LatticeBoltzmann Yagi(Params,Paras,JSources);
  int N=100;
  vector<double> S(N,0),S_current(N,0),Splus(N,0);
  vector<double> S2(N,0),S2_current(N,0),S2plus(N,0);
  double R = 2*lamd;
  int t, tmax=3*T;
  Yagi.Start();

  int Ntheta = 100, Nphi = 100;
  vector<vector<double>> Ss,SsCurrent,SsPlus;
  Set(Ss,Ntheta,Nphi);
  Set(SsCurrent,Ntheta,Nphi);
  Set(SsPlus,Ntheta,Nphi);

  for(t=0;t<tmax;t++){
    Yagi.UpdateTime(t);
    cout<< "Iteracion:\t"<< t+1 <<endl;
    Yagi.Collision();
    Yagi.Advection();

    if(R <= C*t){

      if(Ly==1 || D3d){
        //Yagi.PowerPlanePhi(R,N,S_current);
        Yagi.PowerPlane(R,N,M_PI/4,S_current);
        transform(S.begin(),S.end(),S_current.begin(),Splus.begin(),plus<double>());
        S.assign(Splus.begin(),Splus.end());
        //compare(S,S_current,N);
      }

      if(Lz==1 || D3d){
        //Yagi.PowerPlane(R,N,M_PI/2,S2_current);
        Yagi.PowerPlaneTheta(R,N,S2_current);
        //Yagi.PowerPlanePhi(R,N,S2_current);
        transform(S2.begin(),S2.end(),S2_current.begin(),S2plus.begin(),plus<double>());
        S2.assign(S2plus.begin(),S2plus.end());
        //compare(S2,S2_current,N);
      }

      if(D3d)
        Yagi.PowerByPlanes(R,Ntheta,Nphi,Ss,SsCurrent,SsPlus);

    }

  }

  //Yagi.Print();
  //Print Power data
  Theo formulas; formulas.Init(Params);
  if(Ly==1 || D3d )
    PrintPower(S,N,R,"EPlane.txt",formulas,1);
  if(Lz==1 || D3d)
    PrintPower(S2,N,R,"BPlane.txt",formulas,0);
  //PrintPower(S,N,R,"EPlane.txt",formulas,1);

  if(D3d)
    Yagi.PrintPowerPlanes(Ss,Ntheta,Nphi,R);

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
      file2 << phi  << " " << Pteo/(Z0*J0*J0) << endl;
      file << M_PI/2-phi<< " " << V[i]/(Z0*J0*J0)*0.5 << endl;
    }
    else{
      Pteo = formulas.Power(M_PI/2,phi,R);
      file2 << phi << " " << Pteo/(Z0*J0*J0) << endl;
      file << phi << " " << V[i]/(Z0*J0*J0)*0.5 << endl;
    }

  }
  file.close();
  file2.close();
}
