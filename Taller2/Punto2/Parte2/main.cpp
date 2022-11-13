#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;


//Indices de refraccion

const int Lx=600;
const int Ly=200;

const int Q=5;
const double W0=1.0/3;



const double C=0.5; // C<0.707 cells/click
const double C2=C*C;



const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//-------Clase LatticeGas--------------

class LatticeBoltzmann{
  private:
    double w[Q];      //Weights
    int Vx[Q],Vy[Q];  //Velocity vectors
    double *f, *fnew; //Distribution Functions

  public:
    LatticeBoltzmann();
    ~LatticeBoltzmann();
    int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
    double rho(int ix,int iy,bool UseNew);
    double Jx(int ix,int iy,bool UseNew);
    double Jy(int ix,int iy,bool UseNew);
    double Ccelda(int ix,int iy);
    double feq(double rho0,double Jx0,double Jy0,int i,int ix,int iy);
    void Start(double rho0,double Jx0,double Jy0);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    void Print(const char * NombreArchivo);

};

LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
}

double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int direccion,n0;
  for(sum=0,direccion=0;direccion<Q;direccion++){
    n0=n(ix,iy,direccion);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int direccion,n0;
  for(sum=0,direccion=0;direccion<Q;direccion++){
    n0=n(ix,iy,direccion);
    if(UseNew) sum+=Vx[direccion]*fnew[n0]; else sum+=Vx[direccion]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int direccion,n0;
  for(sum=0,direccion=0;direccion<Q;direccion++){
    n0=n(ix,iy,direccion);
    if(UseNew) sum+=Vy[direccion]*fnew[n0]; else sum+=Vy[direccion]*f[n0];
  }
  return sum;
}

double LatticeBoltzmann::Ccelda(int ix,int iy){
  return C;
}


double  LatticeBoltzmann::feq(double rho0,double Jx0,
                              double Jy0,int i, int ix, int iy){
  double v2 = pow(Ccelda(ix,iy),2);
  double Aux = 1-3*v2*(1-W0);
  if(i>0)
    return 3*w[i]*(v2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return rho0*Aux;
}

void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
        n0=n(ix,iy,i);
        f[n0]=feq(rho0,Jx0,Jy0,i,ix,iy);
      }
}


void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Jx0,Jy0;
  double check_mirror;
  double mx0 = 50, my0 = 100, R = 100; //mirror coordinates
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){

      check_mirror = pow(ix-mx0,2)+pow(iy-my0,2) - R*R;
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++){ //for each velocity vector
        n0=n(ix,iy,i);

        if(ix > 100 && ix <200 && check_mirror>0){
          if(i==0)
            fnew[n0] = f[n(ix,iy,0)];
          else
            fnew[n0] = f[n(ix,iy,(Q-1+i-1+(Q-1)/2)%(Q-1)+1)];
        }else
          fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i,ix,iy);
      }
    }
}

void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0;
  double lambda,omega,rho0,Jx0,Jy0,v;
  //an oscillating source in ix=0
  for(ix=0,iy=0;iy < Ly; iy++){
    v = Ccelda(ix,iy);
    lambda=10; omega=2*M_PI/lambda*v;
    rho0=10*sin(omega*t); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(rho0,Jx0,Jy0,i,ix,iy);
    }
  }
}

void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
        n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
        f[n0next]=fnew[n0]; //periodic boundaries
      }
}

void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix,iy;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//------------Programa Principal-------------

int main(){
  LatticeBoltzmann Ondas;
  int t, tmax = 500;
  double rho0=0,Jx0=0,Jy0=0;
  Ondas.Start(rho0,Jx0,Jy0);

  for(t=0;t<tmax;t++){
    Ondas.Collision();
    Ondas.ImposeFields(t);
    Ondas.Advection();
  }

  //Print
  Ondas.Print("data.dat");

  return 0;
}
