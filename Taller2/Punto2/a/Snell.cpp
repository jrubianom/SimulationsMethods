#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

const int Lx = 400;
const int Ly = 200;

const int ix0 = 100; 
double ThetaInterface_i = 30.;
double ThetaInterface = (90.-ThetaInterface_i)*M_PI/180.;

const int Q =5;
const double W0 = 1./3.;

const double C = 0.5;// C<0.707 cells/click
const double C2 = C*C;
const double AUX0 = 1-3*C2*(1-W0);

const double tau = 0.5;
const double Utau = 1./tau;
const double UmUtau = 1- Utau;

//--------------------Clase Lattice Boltzmann-------------
class LatticeBoltzmann{
private:
   double w[Q];
   double Vx[Q], Vy[Q];
   double *f, *fnew;

public:
   LatticeBoltzmann(void);
   ~LatticeBoltzmann(void);
   int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
   double Ccelda(int ix, int iy);
   double C2celda(int ix, int iy);
   double AUXcelda(int ix, int iy);
   double rho(int ix,int iy,bool UseNew);
   double Jx(int ix,int iy,bool UseNew);
   double Jy(int ix,int iy,bool UseNew);
   double feq(int ix, int iy, double rho0,double Jx0,double Jy0,int i);
   void Start(double rho0,double Jx0,double Jy0);
   void Collision(void);
   void ImposeFields(int t);
   void Advection(void);
   void Print(string NameFile, int Rangex, int Rangey);
   void Perfil(int iy,int Startx,int Rangex);
};

LatticeBoltzmann::LatticeBoltzmann(void){
   //set weights
   w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.-W0)/4.;
   //Set velocity vectors:
   Vx[0]=0; Vx[1]=1; Vx[2]=0; Vx[3]=-1; Vx[4]=0;
   Vy[0]=0; Vy[1]=0; Vy[2]=1; Vy[3]=0; Vy[4]=-1;
   //Create Dynamic arrays
   int ArraySize = Lx*Ly*Q;
   f = new double [ArraySize]; fnew = new double [ArraySize];
}

LatticeBoltzmann::~LatticeBoltzmann(void){
   delete [] f; delete [] fnew;
}

double LatticeBoltzmann::Ccelda(int ix, int iy){
  int ixInterface, iyInterface = Ly/2.0;
  ixInterface = (iy-iyInterface)/(tan(ThetaInterface))+ix0;
  double C = tanh(ixInterface-ix)*(0.125)+(0.125*3);
  return C;
}

double LatticeBoltzmann::C2celda(int ix, int iy){
  double C = Ccelda(ix,iy);
  return C*C;
}
double LatticeBoltzmann::AUXcelda(int ix, int iy){
  double C2 = C2celda(ix,iy);
  double AUX0 = 1-3*C2*(1-W0);
  return AUX0;
}



double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}  
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}  


double  LatticeBoltzmann::feq(int ix, int iy,double rho0,double Jx0,double Jy0,int i){
  if(i>0)
    return 3*w[i]*(C2celda(ix,iy)*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return rho0*AUXcelda(ix,iy);
}  

void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	n0=n(ix,iy,i);
	f[n0]=feq(ix,iy,rho0,Jx0,Jy0,i);
      }
}  
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++){ //for each velocity vector
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*feq(ix,iy,rho0,Jx0,Jy0,i);
      }
    }  
}

void LatticeBoltzmann::ImposeFields(int t){
  int i,ix,iy,n0;
  double lambda,omega,rho0,Jx0,Jy0; lambda=10; omega=2*M_PI/lambda*C;
  //an oscillating source in the middle
  ix = 0;
  rho0=10*sin(omega*t);
  for (iy = 0; iy< Ly; iy++){
   Jx0=Jx(ix,iy,false);
   Jy0=Jy(ix,iy,false);
   for(i=0;i<Q;i++){
     n0=n(ix,iy,i);
     fnew[n0]=feq(ix,iy,rho0,Jx0,Jy0,i);
   }
  }
}
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext=int(ix+Vx[i]+Lx)%int(Lx); iynext=int(iy+Vy[i]+Ly)%int(Ly);
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}
void LatticeBoltzmann::Print(string NameFile, int Rangex, int Rangey){
  ofstream MyFile(NameFile); double rho0; int ix,iy;
  for(ix=0;ix<Rangex;ix++){
    for(iy=0;iy<Rangey;iy++){
      rho0=rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

void LatticeBoltzmann::Perfil(int iy, int Startx, int Rangex){
  int ix;
  
  ofstream perfil ("Perfil"+to_string(iy)+".dat");
  for (ix = Startx; ix< Rangex; ix++){
    perfil <<ix<<"\t"<<rho(ix,iy,true)<<endl;
  }  
  perfil.close();
}
//----------- Funciones Globales -----------

void InicieCuadro(ofstream &file,string frame){
    file<<"splot '" << frame <<"'"<<endl;
}

void TermineCuadro(ofstream &file){
    file<<endl;
}

void InicieAnimacion(ofstream &file, int Rangex, int Rangey){
  file<<"set terminal gif animate"<<endl;
  file<<"set output 'snell.gif'"<<endl;
  file<<"unset key"<<endl;
  file<<"set xrange[0:"<<Rangex<<"]"<<endl;
  file<<"set yrange[0:"<<Rangey<<"]"<<endl;
  file<<"set size ratio -1"<<endl;
  file<<"set pm3d map"<<endl;
  file<<"set cbrange [-25:25]"<<endl;
}


//------------Programa Principal-------------

int main(void){
   LatticeBoltzmann Ondas;
   int t , tmax= 400;
   double rho0 = 0, Jx0 =0, Jy0 =0; 
   
   string current_frame = "Animacion/frame.dat";
   ofstream AnimFile("animation.gp");
   InicieAnimacion(AnimFile,200,200);
   
   Ondas.Start(rho0,Jx0,Jy0);
   
   for (t=0; t<tmax; t++){
      Ondas.Collision();
      Ondas.ImposeFields(t);
      Ondas.Advection();   
      if(t%1 == 0 ){
        Ondas.Print(current_frame,200,200);
        InicieCuadro(AnimFile,current_frame);
        current_frame = "Animacion/frame"+to_string(t)+".dat";
      //TermineCuadro(AnimFile);
      }
   }
   
   Ondas.Print("Interfaz.dat",200,200);
   Ondas.Perfil(100,80,200);
   Ondas.Perfil(0,30,200);
   
   AnimFile.close();
  return 0;
}
