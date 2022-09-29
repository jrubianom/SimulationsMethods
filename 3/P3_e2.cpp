#include <iostream>
#include <fstream>
#include <cmath>
#include "vector.h"

/*
Misma plantilla del inciso D.
El único cambio es que se imprime el tiempo
y la posición x del troyano.
*/

const int N=3;

const double G = 1.0;

//----------Constantes PEFRL
const double zeta = 0.1786178958448091e00;
const double lambda = -0.2123418310626054e00;
const double chi = -0.6626458266981849e-01;

const double coef1 = (1 - 2*lambda)/2;
const double coef2 = 1 - 2*(chi + zeta);

//-----------Declaración de clases
class Cuerpo;
class Interaccionador;

//-----------Clase Cuerpo
class Cuerpo{
private:
  vector3D r, V, F;
  double m, R;
public:
  void inicie(double x0, double y0, double z0,
	      double Vx0, double Vy0, double Vz0,
	      double m0, double R0);
  void borrar_fuerza(void);
  void sumar_fuerza(vector3D F0);
  void mover_r(double dt, double coef);
  void mover_v(double dt, double coef);
  double get_x(void){return r.x();};
  double get_y(void){return r.y();};
  double get_z(void){return r.z();};
  double get_r(void){return r.norm();};
  double get_theta(void);
  friend class Interaccionador;
};

vector3D rotar_marco(Cuerpo C, double theta);

//-----------Clase Interaccionador
class Interaccionador{
public:
  void calcule_fuerzas(Cuerpo *planeta_arr);
  void calcule_fuerzas_entre(Cuerpo &planeta1, Cuerpo &planeta2);
};

//-----------Funciones globales
// Imprime los comandos de gnuplot para producir la imagen
void gnuplot_output();
  
int main() {
  Cuerpo planeta_arr[N];
  Interaccionador Newton;

  double m0=1047, m1=1, r=1000;
  double M=m0+m1, x0=-m1*r/M, x1=m0*r/M;
  double omega=std::sqrt(G*M/(r*r*r)), T=2*M_PI/omega, V0=omega*x0, V1=omega*x1;
  double t, dt=0.01;
  double tdibujo, tcuadro=T/25;
  double tmax=20*T;

  int i;
  
  // Condición inicial del tercer cuerpo
  double m2=0.005;
  double x2=r*std::cos(M_PI/3), y2=r*std::sin(M_PI/3);
  double V2x=-omega*y2, V2y=omega*x2;

  vector3D posicion_rotada[N];
  
  planeta_arr[0].inicie(x0, 0, 0, 0, V0, 0, m0, 1.0);
  
  planeta_arr[1].inicie(x1, 0, 0, 0, V1, 0, m1, 0.5);

  planeta_arr[2].inicie(x2, y2, 0, V2x, V2y, 0, m2, 0.1);

  std::ofstream fout("P3_e.dat");
  fout << t << " " << planeta_arr[2].get_x() << "\n";
  
  for(t=0, tdibujo=0 ; t<tmax; t+=dt, tdibujo+=dt){
    // Rotar marco de referencia
    for(i=0;i<N;i++){
      double theta = planeta_arr[1].get_theta(); // Ángulo de Júpiter
      posicion_rotada[i] = rotar_marco(planeta_arr[i], theta);
    }
    // Guardar datos
    if(tdibujo>tcuadro){
      fout << t << " " << posicion_rotada[2].x() << "\n";
      tdibujo=0;
    }
    
    for(i=0;i<N;i++) planeta_arr[i].mover_r(dt,zeta);
    
    Newton.calcule_fuerzas(planeta_arr);
    for(i=0;i<N;i++) planeta_arr[i].mover_v(dt,coef1);
    
    for(i=0;i<N;i++) planeta_arr[i].mover_r(dt,chi);
    
    Newton.calcule_fuerzas(planeta_arr);
    for(i=0;i<N;i++) planeta_arr[i].mover_v(dt,lambda);

    for(i=0;i<N;i++) planeta_arr[i].mover_r(dt,coef2);

    Newton.calcule_fuerzas(planeta_arr);
    for(i=0;i<N;i++) planeta_arr[i].mover_v(dt,lambda);

    for(i=0;i<N;i++) planeta_arr[i].mover_r(dt,chi);

    Newton.calcule_fuerzas(planeta_arr);
    for(i=0;i<N;i++) planeta_arr[i].mover_v(dt,coef1);

    for(i=0;i<N;i++) planeta_arr[i].mover_r(dt,zeta);
  }
  for(i=0;i<N;i++){
    double theta = planeta_arr[1].get_theta();
    posicion_rotada[i] = rotar_marco(planeta_arr[i], theta);
  }
  fout.close();
  gnuplot_output();
  return 0;
}

//----------Cuerpo
void Cuerpo::inicie(double x0, double y0, double z0, double Vx0, double Vy0, double Vz0, double m0, double R0){
  r.load(x0, y0, z0);
  V.load(Vx0, Vy0, Vz0);
  m = m0; R = R0;
}

void Cuerpo::borrar_fuerza(void){
  F *= 0;
}

void Cuerpo::sumar_fuerza(vector3D F0){
  F += F0;
}

void Cuerpo::mover_r(double dt, double coef){
  r += V*(dt*coef);
}

void Cuerpo::mover_v(double dt, double coef){
  V += (F/m)*(dt*coef);
}

double Cuerpo::get_theta(void){
  return std::atan2(r.y(), r.x());
}

vector3D rotar_marco(Cuerpo C, double theta){
  double x=C.get_x(), y=C.get_y();
  double x_new = std::cos(theta)*x + std::sin(theta)*y;
  double y_new = -std::sin(theta)*x + std::cos(theta)*y;

  vector3D r_new; r_new.load(x_new, y_new, 0);

  return r_new;
}

//----------Interaccionador
void Interaccionador::calcule_fuerzas(Cuerpo *planeta_arr){
  int i,j;
  //Borrar fuerzas
  for(i=0;i<N;i++)
    planeta_arr[i].borrar_fuerza();
  //Calcular fuerzas entre parejas
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++)
      calcule_fuerzas_entre(planeta_arr[i], planeta_arr[j]);
  }
}

void Interaccionador::calcule_fuerzas_entre(Cuerpo &planeta1, Cuerpo &planeta2){
  vector3D r21, F;
  double m2 = planeta2.m, m1 = planeta1.m;

  r21 = planeta2.r - planeta1.r;
  double d = r21.norm();
  
  F = -G*m1*m2*r21/std::pow(d,3);
  planeta2.sumar_fuerza(F);
  planeta1.sumar_fuerza(F*(-1));
}

void gnuplot_output(){
  std::cout << "set terminal pdf" << std::endl;
  std::cout << "set output 'P3_e.pdf'" << std::endl;
  //std::cout << "set size ratio -1" << std::endl;
  std::cout << "plot \"P3_e.dat\" u 1:2 w l t 'Posicion X'" << std::endl;
}
