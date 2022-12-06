#include<iostream>
#include<cmath>

using namespace std;
class Parameter{
    private:
        //------//Dimensions
        int Lx,Ly,Lz;
        //----------- taus
        double Tau = 0.5;
        double UTau = 1/Tau;
        double UmUTau=1-1/Tau;
        //------------
        double Epsilon0=1, Mu0=2;
        double C=1.0/sqrt(Epsilon0*Mu0);
        double Sigma0=0.0,Z0;
        //-----------------
        double E00=0.001,B00=E00/C;
        double J0=0.0001;
        double alpha = 0.25; //alpha = parameter related with width of the Gaussian

        double T=25;
        double k,lambda,omega;

    public:
        Parameter(int Lx0,int Ly0,int Lz0,double Tau0,double Epsilon00,
                  double Mu00,double Sigma00,double Ei,double Bi,
                  double Ji,double alpha0,double T0);
        friend class LatticeBoltzmann;
        friend class Theo;
};

Parameter::Parameter(int Lx0,int Ly0,int Lz0,double Tau0,double Epsilon00,
                     double Mu00,double Sigma00,double Ei,double Bi,
                     double Ji,double alpha0,double T0){
    Lx=Lx0; Ly=Ly0; Lz=Lz0;
    Tau=Tau0; UTau = 1/Tau; UmUTau=1-1/Tau;
    Epsilon0=Epsilon00; Mu0=Mu00;  C=1.0/sqrt(Epsilon0*Mu0); Sigma0=Sigma00;
    E00=Ei; B00=Bi; J0=Ji;
    alpha=alpha0;
    T=T0;
    Z0 = sqrt(Mu0/Epsilon0);
    omega = 2*M_PI/T; lambda = C*T; k = omega/C;
}
