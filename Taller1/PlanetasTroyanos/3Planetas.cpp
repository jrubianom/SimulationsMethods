#include <iostream>
#include <cmath>
#include "vector.h"
#include<random>

using namespace std;

//Constantes globales
const double G = 1.0;
const int N = 3; //Numero de Planetas

//Constamtes de  PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaracion de Clases
class Cuerpo;
class Colisionador;

//-------------------Clase Cuerpo --------------
class Cuerpo{
  private:
    vector3D r,V,F;
    double m,R;
  public:
    void Inicie(double x0,double y0,double z0,
                double Vx0,double Vy0,double Vz0,
                double m0,double R0);
    void BorreFuerza(void) {F.load(0,0,0);};
    void SumeFuerza(vector3D F0){F += F0;};
    void Mueva_r(double dt, double coeff);
    void Mueva_v(double dt, double coeff);
    double Getx(void){return r.x();};
    double Gety(void){return r.y();};
    double Getz(void){return r.z();};
    friend class Colisionador;
    friend vector3D Rotate(double theta,Cuerpo P,char atributo);
    friend double GetTheta(Cuerpo P);
};

//--------------Clase Colisionador------------------
class Colisionador{
  private:
  public:
    void CalculeFuerzas(Cuerpo * Planeta);
    void CalculeFuerzasEntre(Cuerpo & Planeta1,Cuerpo &Planeta2);
};

//-----------------Definicion de metodos de Cuerpo------------
void Cuerpo::Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0,double Vz0,
                    double m0,double R0){
  r.load(x0,y0,z0);
  V.load(Vx0,Vy0,Vz0);
  m = m0; R = R0;
}

void Cuerpo::Mueva_r(double dt, double coeff){
  r += V*(coeff*dt);
}

void Cuerpo::Mueva_v(double dt, double coeff){
  V += F*(dt*coeff/m);
}


void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
  int i,j;
  //Borrar Fuerzas
  for(i=0;i<N;i++){
    Planeta[i].BorreFuerza();
  }
  //Calcular las fuerzas entre todas las parejas de planetas
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      CalculeFuerzasEntre(Planeta[i],Planeta[j]);
    }
  }
}

//Given a vector u in I frame, gets the coordinates of u in the rotated frame
vector3D Rotate(double theta,Cuerpo P,char atributo){
  vector3D R1,R2,rnew;
  R1.load(cos(theta),sin(theta),0);
  R2.load(-sin(theta),cos(theta),0);
  if(atributo == 'r') rnew.load(P.r*R1,P.r*R2,0);
  else if(atributo == 'v') rnew.load(P.V*R1,P.V*R2,0);
  return rnew;
}

double GetTheta(Cuerpo P){
  vector3D e1; e1.load(1,0,0);
  double theta = angle(e1,P.r);
  if(P.Gety() < 0 ) theta *= -1;
  return theta;
}

//--------------------Definicion de metodos de Colisionador--------------
void Colisionador::CalculeFuerzasEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
  vector3D r21,n,F1; double d,F;
  r21 = Planeta2.r - Planeta1.r; d = r21.norm(); n = r21/d;
  F = G*Planeta1.m*Planeta2.m*pow(d,-2.0);
  F1 = F*n;
  Planeta1.SumeFuerza(F1); Planeta2.SumeFuerza(F1*(-1));
}

vector3D random_vector(double norma,double eps = 1e-3){
  vector3D rv;
  const int SEED = 1;
  std::mt19937 gen(SEED); // declarando el generador
  std::uniform_real_distribution<double> dis(-1.0, 1.0);
  rv.load(dis(gen),dis(gen),0); rv /= rv.norm();
  return rv*(eps*norma);
}


int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double m0=1047,m1=1,m2 = 0.005,r=1000;
  double M=m0+m1,x0 = -m1*r/M,x1=m0*r/M;
  double omega = sqrt(G*M/pow(r,3)), T = 2*M_PI/omega, tmax = 30*T;
  double V0 = omega*x0,V1 = omega*x1;
  double t, dt = 0.1;
  double theta = M_PI/3;
  //--------------- Inicializacion Sol y Jupyter (Planetas 0 y 1) ------------
  Planeta[0].Inicie(x0,0,0, 0,V0, 0,m0,1.0);
  Planeta[1].Inicie(x1,0,0, 0,V1,0,m1,0.5);
  //-----------------Inicializacion Trojano (planeta 2) ------------
  vector3D Tr0 = Rotate(-theta,Planeta[1],'r');
  vector3D Tv0 = Rotate(-theta,Planeta[1],'v');
  //vector de perturbacion de velocidad
  vector3D ee = random_vector(V1,0);
  Planeta[2].Inicie(Tr0.x(),Tr0.y(),Tr0.z(),
                    Tv0.x()+ee.x(),Tv0.y()+ee.y(),Tv0.z()+ee.z(),
                    m2,0.005);
  vector3D Tr_non; //Coordenadas de Trojano en el marco no incerical
  int index_planeta = 2;
  //Evolucion dinamica
  int j = 0;
  int step = (int) (T*0.1/dt);
  for(t = 0; t < tmax; t += dt){
    theta = GetTheta(Planeta[1]);
    //theta = 0;
    Tr_non = Rotate(theta,Planeta[index_planeta],'r');
    if(j % step == 0) cout << t/T << "\t" << Tr_non.x() << "\t" << Tr_non.y() << endl;
    // Mover por PEFRL
    for(int i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zeta);
    Newton.CalculeFuerzas(Planeta);
    for(int i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);
    for(int i=0;i<N;i++) Planeta[i].Mueva_r(dt,Chi);
    Newton.CalculeFuerzas(Planeta);
    for(int i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);
    for(int i=0;i<N;i++) Planeta[i].Mueva_r(dt,Coeficiente2);
    Newton.CalculeFuerzas(Planeta);
    for(int i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);
    for(int i=0;i<N;i++) Planeta[i].Mueva_r(dt,Chi);
    Newton.CalculeFuerzas(Planeta);
    for(int i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);
    for(int i=0;i<N;i++) Planeta[i].Mueva_r(dt,Zeta);
    j++;
  }

  return 0;
}
