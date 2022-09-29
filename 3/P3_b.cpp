#include <iostream>
#include <fstream>
#include <cmath>
#include "vector.h"

/*
La integración del sistema se realiza de la misma forma
que en inciso A. Para obtener la posición en el marco
rotado, tras cada paso de PEFRL se calcula el ángulo que
el Sol y Júpiter se traslada. Luego, los vectores de
posición se rotan la cantidad correspondiente. El resultado
de la rotación se almacena en el arreglo posicion_rotada.
Para la rotación se crea la función rotar_marco. Note que
esta no modifica los vectores del argumento.
*/

const int N=2;

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

vector3D rotar_marco(Cuerpo C);

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
  double tdibujo, tcuadro=T/100;
  double tmax=20*T;
  
  int i;

  vector3D posicion_rotada[N];
  
  planeta_arr[0].inicie(x0, 0, 0, 0, V0, 0, m0, 1.0);
  
  planeta_arr[1].inicie(x1, 0, 0, 0, V1, 0, m1, 0.5);

  std::ofstream fout("P3_b.dat");
  fout << planeta_arr[0].get_x() << " " << planeta_arr[0].get_y() << " "
       << planeta_arr[1].get_x() << " " << planeta_arr[1].get_y() << "\n";
  
    
  for(t=0, tdibujo=0 ; t<tmax; t+=dt, tdibujo+=dt){
    // Rotar marco de referencia
    for(i=0;i<N;i++) posicion_rotada[i] = rotar_marco(planeta_arr[i]);
    
    // Guardar datos
    if(tdibujo>tcuadro){
      fout << posicion_rotada[0].x() << " " << posicion_rotada[0].y() << " "
	   << posicion_rotada[1].x() << " " << posicion_rotada[1].y() << "\n";
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
  for(i=0;i<N;i++) posicion_rotada[i] = rotar_marco(planeta_arr[i]);

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

vector3D rotar_marco(Cuerpo C){
  double x=C.get_x(), y=C.get_y(), theta=C.get_theta();
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
  std::cout << "set output 'P3_b.pdf'" << std::endl;
  //std::cout << "set size ratio -1" << std::endl;
  std::cout << "plot \"P3_b.dat\" u 1:2 w l, \"\" u 3:4 w l" << std::endl;
}
