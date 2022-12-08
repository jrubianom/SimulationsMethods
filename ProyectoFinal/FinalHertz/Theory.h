#include<iostream>
#include<cmath>
#include "Parameters.h"

using namespace std;
class Theo{
        private:
                int Lx,Ly,Lz;
                //----------- taus
                double Tau, UTau,UmUTau;
                //------------
                double Epsilon0,Mu0;
                double C,Sigma0,Z0;
                //-----------------
                double E00,B00,J0;
                double alpha; //alpha = parameter related with width of the Gaussian
                double T,lambda,k,omega;
                double p;
        public:
                void Init(Parameter P0);
                friend class LatticeBoltzmann;
                double B_y(double r,double t);
                double E_z(double r,double t);
                //double AmplitudeB(double r,double t=0);
                //double AmplideE(double r,double t=0);
                double Power(double theta,double phi,double R);
};

void Theo::Init(Parameter P0){
        Lx = P0.Lx; Ly = P0.Ly; Lz = P0.Lz;
        Tau = P0.Tau; UTau = P0.UTau; UmUTau = P0.UmUTau;
        Epsilon0 = P0.Epsilon0; Mu0 = P0.Mu0; C=P0.C; Sigma0 = P0.Sigma0;
        Z0 = sqrt(Mu0/Epsilon0);
        E00=P0.E00; B00=P0.B00;
        J0=P0.J0;
        alpha=P0.alpha;
        T = P0.T;
        omega = 2*M_PI/T; lambda = C*T; k = omega/C;
        p = J0/omega*pow(M_PI/alpha,1.5);
}


double Theo::B_y(double r,double t){
        double A = k*k*Z0*p/(4*M_PI);
        return A/r*(cos(k*r-omega*t)-sin(k*r-omega*t)/(k*r));
}

double Theo::E_z(double r,double t){
        double A = C*k*k*Z0*p/(4*M_PI);
        return A/r*(-sin(k*r-omega*t)/(k*r) + cos(k*r-omega*t)*(1-1.0/pow(k*r,2)));
}

//Get max E X B
double Theo::Power(double theta,double phi,double R){
        double Ab = k*k*Z0*p/(4*M_PI);
        double Ae = C*Ab;
        return Ab*Ae*pow(sin(theta),2)*0.5;
}
