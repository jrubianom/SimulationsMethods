#include <iostream>
#include <cmath>
#include "../../Vector3D/vector.h"

//Constantes globales
const int N=2;
const double G=1.0;

//consantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Clases
class Cuerpo;
class Colisionador;

class Cuerpo{
private:
	vector3D r,V,F;
	double m,R;
public:
	void Init(double x0,double y0, double z0,
			  double Vx0,double Vy0, double Vz0,
			  double m0,double R0);
	void Force_erase();
	void Force_add(vector3D F0);
	void Move_r(double dt, double coeff);
	void Move_V(double dt, double coeff);
	double Getx(){return r.x();}
	double Gety(){return r.y();}
	double Getz(){return r.z();}
	double Get_rnorm(){return r.norm();}
	void Dibujese(void);
	friend class Colisionador;
};
void Cuerpo::Init(double x0,double y0, double z0,double Vx0,double Vy0, double Vz0,double m0,double R0){
	r.load(x0,y0,z0); V.load(Vx0,Vy0,Vz0);
	m=m0; R=R0;
}

void Cuerpo::Force_erase(){
	F.load(0,0,0);
}

void Cuerpo::Force_add(vector3D F0){
	F+=F0;
}

void Cuerpo::Move_r(double dt, double coeff){
	r+=V*(dt*coeff);
}

void Cuerpo::Move_V(double dt, double coeff){
	V+=F*(dt*coeff/m);
}

void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

class Colisionador{
private:
public:
	void Force(Cuerpo *Planeta);
	void Force_betw(Cuerpo & Planeta1, Cuerpo & Planeta2);
};

void Colisionador::Force(Cuerpo *Planeta){
	for(int i=0;i<N;i++) Planeta[i].Force_erase();
	for(int i=0;i<N;i++){
		for(int j=i+1;j<N;j++){
			Force_betw(Planeta[i], Planeta[j]);
		}
	}
}

void Colisionador::Force_betw(Cuerpo & Planeta1, Cuerpo & Planeta2){
	vector3D r21,n,F1; double d21,F;
	r21=Planeta2.r-Planeta1.r; d21=r21.norm(); n=r21/d21;
	F=G*Planeta1.m*Planeta2.m*pow(d21, -2.);
	F1=F*n; Planeta1.Force_add(F1); Planeta2.Force_add(F1*(-1));
}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  //  std::cout<<"set terminal gif animate"<<std::endl;
  //  std::cout<<"set output 'DosPlanetas.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set xrange[-1000:1000]"<<std::endl;
  std::cout<<"set yrange[-1000:1000]"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange [0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
    std::cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    std::cout<<std::endl;
}

// void Graficar_rotado(Cuerpo *planetas){
// 	double x_jupiter,y_jupiter,xpos,ypos;
// 	x_jupiter=planeta[1].Getx();
// 	y_jupiter=(planeta[1]
// 	for(int i=0; i<N; i++){
// 		xpos=(planeta)
// 	}
// }


int main(){
	Cuerpo planeta[N];
	Colisionador Newton;
	double m0=1, m1=1047, r=1000;
	double M=m0+m1, x0=m1*r/M, x1=-m0*r/M;
	double w=sqrt(G*M*pow(r,-3)), T=2*M_PI/w, V0=w*x0, V1=w*x1;
	double t,tmax=20*T, dt=5;
	double tdibujo,tcuadro=T/5000;
	int i;

	planeta[0].Init(x0,0,0,0,V0,0,m0,50);
	planeta[1].Init(x1,0,0,0,V1,0,m1,25);

	InicieAnimacion();
	for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
		//Dibujar animación
		if(tdibujo>tcuadro){
			InicieCuadro();
			for(i=0;i<N;i++) planeta[i].Dibujese();
			TermineCuadro();
			tdibujo=0;
		}

		//Dibujar trayectoria
		// std::cout << planeta[0].Getx() <<"\t"<< planeta[0].Gety() <<"\t"
		// 		  << planeta[1].Getx() <<"\t"<< planeta[1].Gety()
		// 		  << std::endl;

		//Move by PEFRL
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Zeta);
		Newton.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Coeficiente1);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Chi);
		Newton.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Lambda);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Coeficiente2);
		Newton.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Lambda);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Chi);
		Newton.Force(planeta);
		for(i=0;i<N;i++) planeta[i].Move_V(dt,Coeficiente1);
		for(i=0;i<N;i++) planeta[i].Move_r(dt,Zeta);
	}
	return 0;
}
