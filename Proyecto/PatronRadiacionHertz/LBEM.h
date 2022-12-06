#include<iostream>
#include<cmath>
#include "Vector.h"
#include <fstream>
#include "Parameters.h"
#include <vector>
#include <algorithm>

using namespace std;

double Interpolate(double q1,double q2,double q3,double q4,
                  double q5,double q6,double q7,double q8,
                   double x,double y,double z);
double Interpolatebi(double q1,double q2,double q3,double q4,
                     double x,double y);

void SphericalCoordinates(double R,double theta,double phi,double &x,double &y, double &z);


//--------------------- class LatticeBoltzmann ------------
class LatticeBoltzmann{
  private:
    int V[3][3][4], V0[3]; /*V[xyz][p][i]*/  vector3D v[3][4],v0; //v[p][i]
    vector3D e[3][4][2], e0; //e[p][i][j]
    vector3D b[3][4][2], b0; //b[p][i][j]
    double *f=nullptr,*fnew=nullptr;//f[ix][iy][iz][r][p][i][j]
    double *f0 = nullptr,*f0new = nullptr;//f0[ix][iy][iz] (r=0)
    int Qr = 2, Qp = 3, Qi = 4, Qj = 2; // numbers of total r,p,i,j
    //------//Dimensions
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

  public:
    LatticeBoltzmann(Parameter Params0);
    ~LatticeBoltzmann();
    double mur(int ix,int iy,int iz){ return 1.0;}
    double epsilonr(int ix,int iy,int iz){return 1.0;}
    double sigma(int ix,int iy,int iz){return Sigma0;}
    int index(int ix,int iy,int iz,int r,int p,int i,int j);
    int index0(int ix,int iy,int iz);
    //Auxiliary variables
    double prefactor(double epsilonr0,double sigma0){return sigma0/(1+(sigma0*Mu0)/(4*epsilonr0));};
    //Fields from direct sums
    double rhoc(int ix,int iy,int iz,bool UseNew);
    vector3D D(int ix,int iy,int iz,bool UseNew);
    vector3D B(int ix,int iy,int iz,bool UseNew);
    //Fields deduced from the first ones through electromagnetic constants
    vector3D E(vector3D & D0,double Epsilonr);
    vector3D H(vector3D & B0,double Mur);
    //Forced (Actual) Macroscopic Fields (inline)
    vector3D Jprima(vector3D & E0,double prefactor0){return E0*prefactor0;};
    vector3D JprimaSource(int ix,int iy,int iz,int t);
    vector3D Eprima(vector3D & E0,vector3D & Jprima0,double epsilonr0){return E0-Jprima0*(Mu0/(4*epsilonr0));};
    vector3D Dprima(vector3D & Eprima0,double epsilonr0){return Eprima0/epsilonr0;};
    //Equilibrium Functions
    double feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
               double Epsilonr,double Mur,
               int r,int p,int i,int j);
    double feq0(double rhoc0);
    //Simulation Functions
    void Start(void);
    void Collision(int t);
    void ImposeFields(int t);
    void Advection(void);
    void MacroscopicFields(int ix,int iy, int iz,
                           vector3D &Eprima0,vector3D &H0, vector3D &B0);
    double Poynting(int ix,int iy,int iz);
    double PowerAtPoint(double x,double y, double z);
    void PowerOverSphere(double R,int Ntheta,int Nphi,vector<double> &V);
    double PowerAtPointTheta(double x,double y);
    double PowerAtPointPhi(double x,double z);
    void PowerPlanePhi(double R,int Ntheta,vector<double> &V);
    void PowerPlaneTheta(double R,int Nphi,vector<double> &V);
    void Print(void);
};


LatticeBoltzmann::LatticeBoltzmann(Parameter Params0){

  Lx = Params0.Lx; Ly = Params0.Ly; Lz = Params0.Lz;
  Tau = Params0.Tau; UTau = Params0.UTau; UmUTau = Params0.UmUTau;
  Epsilon0 = Params0.Epsilon0; Mu0 = Params0.Mu0; C=Params0.C; Sigma0 = Params0.Sigma0;
  Z0 = sqrt(Mu0/Epsilon0);
  E00=Params0.E00; B00=Params0.B00;
  J0=Params0.J0;
  alpha=Params0.alpha;
  T = Params0.T;
  omega = 2*M_PI/T; lambda = C*T; k = omega/C;

  int ix,iy,iz,alpha,r,p,i,j;
  //Velocity vectors V[p][i]=V^p_i (in components)
  V0[0]=V0[1]=V0[2]=0;

  V[0][0][0]=V[0][1][0]=V[1][2][0]=1;
  V[1][0][0]=V[2][1][0]=V[2][2][0]=1;
  V[2][0][0]=V[1][1][0]=V[0][2][0]=0;

  V[0][0][1]=V[0][1][1]=V[1][2][1]=-1;
  V[1][0][1]=V[2][1][1]=V[2][2][1]=1;
  V[2][0][1]=V[1][1][1]=V[0][2][1]=0;

  V[0][0][2]=V[0][1][2]=V[1][2][2]=-1;
  V[1][0][2]=V[2][1][2]=V[2][2][2]=-1;
  V[2][0][2]=V[1][1][2]=V[0][2][2]=0;

  V[0][0][3]=V[0][1][3]=V[1][2][3]=1;
  V[1][0][3]=V[2][1][3]=V[2][2][3]=-1;
  V[2][0][3]=V[1][1][3]=V[0][2][3]=0;
  //Velocity vectors V[p][i]=V^p_i (as vectors)
  v0.cargue(V0[0],V0[1],V0[2]); //cargue= load (in Spanish)
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      v[p][i].cargue(V[0][p][i],V[1][p][i],V[2][p][i]);
    }
  //Electric vectors e[p][i][j]=e^p_{ij}
  e0.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++){
      e[p][i][0]=v[p][(i+1)%4]*0.5;
      e[p][i][1]=v[p][(i+3)%4]*0.5;
    }
  //Magnetic vectors b[p][i][j]=b^p_{ij}=v^p_i x e^p_{ij}
  b0.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++)
        b[p][i][j]=(v[p][i]^e[p][i][j]);

  f = new double[Lx*Ly*Lz*Qr*Qp*Qi*Qj]; fnew = new double[Lx*Ly*Lz*Qr*Qp*Qi*Qj];
  f0 = new double[Lx*Ly*Lz]; f0new=new double[Lx*Ly*Lz];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
  delete[] f0; delete[] f0new;
}

int LatticeBoltzmann::index(int ix,int iy,int iz,int r,int p,int i,int j){
  return (iz*Lx*Ly+iy*Lx+ix)*Qr*Qp*Qi*Qj + (r*Qp*Qi*Qj + p*Qi*Qj + i*Qj +j);
}

int LatticeBoltzmann::index0(int ix,int iy,int iz){
  return (iz*Lx*Ly+iy*Lx+ix);
}

//-----------------MACROSCOPIC FIELDS------------------
//Fields from direct sums
double LatticeBoltzmann::rhoc(int ix,int iy,int iz,bool UseNew){
  int p,i,j; double sum;
  int id0,id;
  id0 = index0(ix,iy,iz);
  //Start for the distribution for the central (zero) vector
  if(UseNew)
    sum=f0new[id0];
  else
    sum=f0[id0];
  //Add all the others
  for(p=0;p<2;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,0,p,i,j);
        if(UseNew)
          sum+=fnew[id];
        else
          sum+=f[id];
      }
  return sum;
}
vector3D LatticeBoltzmann::D(int ix,int iy,int iz,bool UseNew){
  int p,i,j,id; vector3D sum; sum.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,0,p,i,j);
        if(UseNew)
          sum+=e[p][i][j]*fnew[id];
        else
          sum+=e[p][i][j]*f[id];
      }
  return sum;
}
vector3D LatticeBoltzmann::B(int ix,int iy,int iz,bool UseNew){
  int p,i,j,id; vector3D sum; sum.cargue(0,0,0);
  for(p=0;p<3;p++)
    for(i=0;i<4;i++)
      for(j=0;j<2;j++){
        id = index(ix,iy,iz,1,p,i,j);
        if(UseNew)
          sum+=b[p][i][j]*fnew[id];
        else
          sum+=b[p][i][j]*f[id];
      }
  return sum;
}
//Fields deduced from the first ones through electromagnetic constants
vector3D LatticeBoltzmann::E(vector3D & D0,double Epsilonr){
  return D0*(1.0/Epsilonr);
}
vector3D LatticeBoltzmann::H(vector3D & B0,double Mur){
  return B0*(1.0/Mur);
}

vector3D LatticeBoltzmann::JprimaSource(int ix,int iy,int iz,int t){
  vector3D Jprima0;
  double Jz;
  Jz=J0*exp(-alpha*(pow((ix-Lx/2),2)+pow((iy-Ly/2),2)+pow((iz-Lz/2),2)))*sin(omega*t);
  Jprima0.cargue(0,0,Jz);
  return Jprima0;
}

//---------------EQUILIBRIUM FUNCTIONS-------------
double LatticeBoltzmann::feq(vector3D & Jprima0,vector3D & Eprima0,vector3D & B0,
                             double epsilonr0,double mur0,
                             int r,int p,int i,int j){
  double VdotJp=(v[p][i]*Jprima0),Epdote=(e[p][i][j]*Eprima0),Bdotb=(b[p][i][j]*B0),aux;
  if(r==0)
    aux=0.25*(0.25*VdotJp+epsilonr0*Epdote+0.5/mur0*Bdotb);
  if(r==1)
    aux=0.25*(0.25*VdotJp+Epdote+0.5*Bdotb);
  return aux;
}
double LatticeBoltzmann::feq0(double rhoc0){
  return rhoc0;
}

//-------------------SIMULATION FUNCTIONS ----------------------------
void LatticeBoltzmann::Start(void){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        prefactor0=prefactor(epsilonr0,sigma0);
        //Impose the fields
        rhoc0=0; D0.cargue(0,0,0); B0.cargue(0,0,0);
        E0=E(D0,epsilonr0); H0=H(B0,mur0);
        Jprima0.cargue(0,0,0);//;=Jprima(E0,prefactor0);
        Eprima0=Eprima(E0,Jprima0,epsilonr0);
        //Impose f=fnew=feq with the desired fields
        id0 = index0(ix,iy,iz);
        f0new[id0]=f0[id0]=feq0(rhoc0);
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                id = index(ix,iy,iz,r,p,i,j);
                fnew[id]=f[id]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              }
      }
}

void LatticeBoltzmann::Collision(int t){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){
        //Compute the constants
        sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
        prefactor0=prefactor(epsilonr0,sigma0);
        //Compute the fields
        rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
        E0=E(D0,epsilonr0); H0=H(B0,mur0);
        //Jprima0=Jprima(E0,prefactor0);
        Jprima0=JprimaSource(ix,iy,iz,t);
        Eprima0=Eprima(E0,Jprima0,epsilonr0);
        //BGK evolution rule
        id0 = index0(ix,iy,iz);
        f0new[id0]=UmUTau*f0[id0]+UTau*feq0(rhoc0);
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                id = index(ix,iy,iz,r,p,i,j);
                fnew[id]=UmUTau*f[id]+UTau*feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
              }
      }
}

/*
  void LatticeBoltzmann::ImposeFields(int t){
  int ix,iy,iz,r,p,i,j; double sigma0,mur0,epsilonr0,prefactor0;
  int id0,id;
  double rhoc0; vector3D D0,B0,E0,H0,Jprima0,Eprima0;
  double T=25.0,omega=2*M_PI/T;//25<T<50
  iz=0; //On the whole incident plane

  for(ix=0;Lx/2-15 < ix && ix<Lx/2 + 15;ix++) //for each cell
  for(iy=0;Ly/2-15 < iy && iy<Ly/2 + 15;iy++)
  for(iz=0;Lz/2-15 < iz && iz<Lz/2 + 15;iz++){
  //Compute the constants
  sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
  prefactor0=prefactor(epsilonr0,sigma0);
  //Compute the fields
  //Primary fields (i.e. those from the sums)
  rhoc0=rhoc(ix,iy,iz,false); D0=D(ix,iy,iz,false); B0=B(ix,iy,iz,false);
  E0=E(D0,epsilonr0); H0=H(B0,mur0);
  Jprima0=JprimaSource(ix,iy,iz,t); Eprima0=Eprima(E0,Jprima0,epsilonr0);
  //Impose fnew=feq with the desired fields
  for(r=0;r<2;r++)
  for(p=0;p<3;p++)
  for(i=0;i<4;i++)
  for(j=0;j<2;j++){
  id0 = index0(ix,iy,iz); id = index(ix,iy,iz,r,p,i,j);
  fnew[id]=feq(Jprima0,Eprima0,B0,epsilonr0,mur0,r,p,i,j);
  f0new[id0]=feq0(rhoc0);
  }
  }
  }
*/

void LatticeBoltzmann::Advection(void){
  int ix,iy,iz,r,p,i,j,ixnew,iynew,iznew;
  int id0,id,idnew;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(iz=0;iz<Lz;iz++){//for each cell
        for(r=0;r<2;r++)
          for(p=0;p<3;p++)
            for(i=0;i<4;i++)
              for(j=0;j<2;j++){
                ixnew=(ix+V[0][p][i]+Lx)%Lx; iynew=(iy+V[1][p][i]+Ly)%Ly; iznew=(iz+V[2][p][i]+Lz)%Lz;
                id0=index0(ix,iy,iz);
                id = index(ix,iy,iz,r,p,i,j);idnew=index(ixnew,iynew,iznew,r,p,i,j);
                f[idnew]=fnew[id];
                f0new[id0]=f0[id0];
              }
      }
}

void LatticeBoltzmann::Print(void){
  int ix,iy=Ly/2,iz,r,p,i,j;
  double rhoc0; vector3D B0,Eprima,H0; double E2,B2;
  ofstream file("Datos/data.txt");
  ofstream file2("Datos/datazeL2.txt");
  for(ix=0;ix < Lx; ix++)
    for(iz=0;iz<Lz;iz++){
      //Compute the electromagnetic constants
      MacroscopicFields(ix,iy,iz,Eprima,H0,B0);
      //Print
      B2 = B0.y()/J0;
      file<<ix<<" "<<iz<<" "<<B2<<endl;
      if(iz==Lz/2)
        file2 << ix << " " << B2 << endl;
    }
  file.close();
  file2.close();
}



//Return the Macroscopic Fields of E,H and B evaluated in a lattice point
void LatticeBoltzmann::MacroscopicFields(int ix,int iy, int iz,
                                         vector3D &Eprima0,vector3D &H0, vector3D &B0){
  double sigma0,mur0,epsilonr0,prefactor0;
  double rhoc0; vector3D D0,E0,Jprima0;

  //Compute the electromagnetic constants
  sigma0=sigma(ix,iy,iz); mur0=mur(ix,iy,iz); epsilonr0=epsilonr(ix,iy,iz);
  prefactor0=prefactor(epsilonr0,sigma0);
  //Compute the Fields
  rhoc0=rhoc(ix,iy,iz,true); D0=D(ix,iy,iz,true); B0=B(ix,iy,iz,true);
  E0=E(D0,epsilonr0); H0=H(B0,mur0);
  Jprima0=Jprima(E0,prefactor0); Eprima0=Eprima(E0,Jprima0,epsilonr0);
}


//Stores in a vector V the values of S*r_hat of a finte set of points over a sphere or radius R
void LatticeBoltzmann::PowerOverSphere(double R,int Ntheta,int Nphi,vector<double> &V){
  double theta,phi,dtheta,dphi;
  double q1,q2,q3,q4,q5,q6,q7,q8;
  double x,y,z;
  dtheta = M_PI/(Ntheta-1); dphi = 2*M_PI/(Nphi-1);
  for(int i=0; i < Nphi; i++){
    phi = i*dphi;
    for(int j=0; j < Ntheta; j++){
      theta=j*dtheta;
      SphericalCoordinates(R,theta,phi,x,y,z);
      x += Lx/2; y += Ly/2; z += Lz/2; //centered in the dipole coordinates
      V[i*Ntheta+j] = PowerAtPoint(x,y,z)/(Z0*J0*J0);
    }
  }
}


void LatticeBoltzmann::PowerPlaneTheta(double R,int Nphi,vector<double> &V){
  double phi, dphi = 2*M_PI/(Nphi-1);
  double x,y,z;
  for(int i=0; i < Nphi; i++){
    phi = i*dphi;
    SphericalCoordinates(R,M_PI/2.0,phi,x,y,z);
    x += Lx/2; y += Ly/2;
    V[i] = PowerAtPointTheta(x,y)/(Z0*J0*J0);
  }
}

void LatticeBoltzmann::PowerPlanePhi(double R,int Ntheta,vector<double> &V){
  double phi, dphi = 2*M_PI/(Ntheta-1);
  double x,y,z;
  for(int i=0; i < Ntheta; i++){
    phi = i*dphi;
    SphericalCoordinates(R,M_PI/2.0,phi,x,z,y);
    x += Lx/2; z += Lz/2;
    V[i] = PowerAtPointPhi(x,z)/(Z0*J0*J0);
  }
}



//Returns the interpolated value of S* r_hat of a point x,y,z
double LatticeBoltzmann::PowerAtPoint(double x,double y, double z){
  int x1 = (int) x, y1 = (int) y, z1 = (int) z;
  double s1,s2,s3,s4,s5,s6,s7,s8;
  s1 = Poynting(x1,y1,z1);
  s2 = Poynting(x1+1,y1,z1);
  s3 = Poynting(x1+1,y1+1,z1);
  s4 = Poynting(x1,y1+1,z1);
  s5 = Poynting(x1,y1,z1+1);
  s6 = Poynting(x1+1,y1,z1+1);
  s7 = Poynting(x1+1,y1+1,z1+1);
  s8 = Poynting(x1,y1+1,z1+1);
  return Interpolate(s1,s2,s3,s4,s5,s6,s7,s8,x,y,z);
}

double LatticeBoltzmann::PowerAtPointTheta(double x,double y){
  int x1 = (int) x, y1 = (int) y;
  double s1,s2,s3,s4,s5,s6,s7,s8;
  s1 = Poynting(x1,y1,Lz/2);
  s2 = Poynting(x1+1,y1,Lz/2);
  s3 = Poynting(x1,y1+1,Lz/2);
  s4 = Poynting(x1+1,y1+1,Lz/2);
  return Interpolatebi(s1,s2,s3,s4,x,y);
}

double LatticeBoltzmann::PowerAtPointPhi(double x,double z){
  int x1 = (int) x, z1 = (int) z;
  double s1,s2,s3,s4,s5,s6,s7,s8;
  s1 = Poynting(x1,Ly/2,z1);
  s2 = Poynting(x1+1,Ly/2,z1);
  s3 = Poynting(x1,Ly/2,z1+1);
  s4 = Poynting(x1+1,Ly/2,z1+1);
  return Interpolatebi(s1,s2,s3,s4,x,z);
}




//--------------Funciones auxuliares
//given a lattice point x,y,z gives the dot product of the poynting vector S* r_hat
double LatticeBoltzmann::Poynting(int ix,int iy,int iz){
  vector3D E,H,B;
  MacroscopicFields(ix,iy,iz,E,H,B);
  vector3D S = E^H;
  return norma(S);
}

//Gets the closest lattice point to (x,y,z)
void closest(double x,double y,double z,int &ix,int &iy, int &iz){
  int x0 = (int) x,y0 = (int) y,z0 = (int) z;
  vector3D r; r.cargue(x,y,z);
  vector3D r0; r0.cargue(x0,y0,z0);
  vector3D diff; diff = r-r0;
  double dist0 = norma(r-r0),dist;
  for(int i=0;i < 2; i++)
    for(int j=0; j < 2; j++)
      for(int k=0; k < 2;k++){
        r0.cargue(x0+i,y0+j,z0+k);
        dist = norma(r-r0);
        if(dist0 > dist){
          ix = x0+i; iy = y0+j; iz = z0+k;
          dist0 = dist;
        }
      }
}

//Makes an trilinear interpolation at point x,y,z.
//q_i is the field evaluated in the i-th point of the cube
double Interpolate(double q1,double q2,double q3,double q4,
                  double q5,double q6,double q7,double q8,
                   double x,double y,double z){
  /*_              8-------------7
    _             -|            -|
    _           -  |           - |
    _         -    |          -  |
    _       5----------------6   |
    _       |     4---------|----3
    _       |    -          |    -
    _       |   -           |   -
    _       |  -            |  -
    _       |-              | -
    _      1----------------2

  */


  int x0 = (int) x, y0 = (int) y, z0 = (int) z;
  double xd = (x-x0), yd = (y-y0), zd = (z-z0);
  double c00,c01,c10,c11;
  c00 = q1*(1-xd)+q2*xd;
  c01 = q5*(1-xd)+q6*xd;
  c10 = q4*(1-xd)+q3*xd;
  c11 = q8*(1-xd)+q7*xd;

  double c0,c1,c;
  c0 = c00*(1-yd)+c10*yd;
  c1 = c01*(1-yd)+c11*yd;
  c = c0*(1-zd)+c1*zd;

  return c;
}

//Transforms from spherical cooridnates to cartesian ones
void SphericalCoordinates(double R,double theta,double phi,double &x,double &y, double &z){
  x = R*sin(theta)*cos(phi);
  y = R*sin(theta)*sin(phi);
  z = R*cos(theta);
}


double Interpolatebi(double q1,double q2,double q3,double q4,
                     double x,double y){
  int ix = (int) x, iy = (int) y;
  double u = x-ix, v = y-iy;
  double qxy = q1*(1-u)*(1-v)+q2*u*(1-v)+q3*(1-u)*v+q4*u*v;
  return qxy;
}
