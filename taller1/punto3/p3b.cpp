#include <iostream>
#include <cmath>
#include "vector.h"

using std::cout;
using std::endl;

//------------------Constantes globales--------------------
const double G = 1.;
const int N = 2,n = 1;

//------------------Constante de PEFRL---------------------
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//------------------Declarar Clases-------------------------

class Cuerpo;
class Colisionador;


//--------------------Clase Cuerpo---------------------------

class Cuerpo{
private: 
	vector3D X, X_Old, V, F;
	double m,R;
public: 
	void BorreFuerza(void){F.load(0,0,0);};
	void SumeFuerza(vector3D F0){F += F0;};
	vector3D GetX(void){return X;};
	vector3D GetV(void){return V;};
	vector3D GetF(void){return F;};
	void Inicie(vector3D X0, vector3D V0, double m0, double R0);
	void Mueva_X(double dt, double coeficiente);
	void Mueva_V(double dt, double coeficiente);
	void Dibujese(void);
	friend class Colisionador;
    	friend double polar(Cuerpo P);
    	friend vector3D rot(Cuerpo C, double theta, char attribute);
};
//----------------funciones de la clase----------------

void Cuerpo::Inicie(vector3D X0, vector3D V0, double m0, double R0){
	X = X0; V = V0;
	m = m0; R = R0;
}

void Cuerpo::Mueva_X(double dt, double coeficiente){
    X += V*(dt*coeficiente);
}

void Cuerpo::Mueva_V(double dt, double coeficiente){
    V += F*(dt*coeficiente)/m;
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<X.x()<<"+"<<R<<"*cos(t),"<<X.y()<<"+"<<R<<"*sin(t)";
}

//---------------------Clase Colisionador--------------------

class Colisionador{
private:

public:
	void CalculeFuerzas(Cuerpo * Planeta);
	void CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2);
};
//--------------funciones colisionador-----------------------

void Colisionador::CalculeFuerzas(Cuerpo * Planeta){
	int i,j;
	//Borrar Fuerzas
	for(i=0;i<N; i++){
		Planeta[i].BorreFuerza();
	}
	//Calcular fuerzas entre todas parejas de planetas
	for(i=0; i<N;i++){
		for(j=i+1; j<N; j++){
			CalculeFuerzaEntre(Planeta[i],Planeta[j]);
		}
	}
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & Planeta1, Cuerpo & Planeta2){
	vector3D r21, n, F1; double d21, F; 
	r21 = Planeta2.X-Planeta1.X; d21 = sqrt(r21.norm2());
	n=r21/d21;  F=G*Planeta1.m*Planeta2.m*pow(d21,-2.0);
  	F1=F*n; Planeta1.SumeFuerza(F1); 
  	Planeta2.SumeFuerza(F1*(-1)); 
}
//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  //cout<<"set terminal gif animate"<<endl; 
  //cout<<"set output 'Jupiter.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-1050:1050]"<<endl;
  cout<<"set yrange[-1050:1050]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}

vector3D rot(Cuerpo C, double theta, int tipo){
	vector3D newX, vec;
	double xprime,yprime;
	if (tipo == 0){ //x
		vec = C.GetX();
	}
	else if (tipo == 1){ //v
		vec = C.GetV();
	}
	xprime = vec.x()*cos(theta)-vec.y()*sin(theta);
	yprime = vec.y()*cos(theta)+vec.x()*sin(theta);
	newX.load(xprime, yprime, vec.z());
	return newX;
}

double polar(Cuerpo P){
  vector3D e1; e1.load(1,0,0);
  double theta = angle(e1,P.X);
  if(P.GetX().y() < 0 ) theta *= -1;
  return theta;
}


//---------------Funciones Globales--------------------------
int main(){
	Cuerpo Planeta[N+n];
  	Colisionador Newton;
  	double m0=1047, m1=1, m2 = 0.005, r=1000;
  	double M=m0+m1, x0=-m1*r/M, x1=m0*r/M;
  	double omega=sqrt(G*M/(r*r*r)), T=2*M_PI/omega, 	v0=omega*x0, v1=omega*x1;
  	double t,tmax=20.1*T,dt=0.1;
  	int i;	
  	
  	double angtroy[3] = {M_PI/3,-M_PI/3, M_PI};
	//-----vectores iniciales----------------
	vector3D X0, X1, V0, V1;
	
	X0.load(x0, 0      ,  0); X1.load(x1, 0     ,  0);
	V0.load(0 , v0,   0); V1.load(0 , v1,  0);
	//----------------------------------------
	
	Planeta[0].Inicie(X0, V0, m0, 50);//sol
  	Planeta[1].Inicie(X1, V1, m1, 55 );//jupiter
  	//------------planeta troyano ------------
  	vector3D Xtroy, Vtroy;
  	
  	Xtroy = rot(Planeta[1], angtroy[0],0);
  	Vtroy = rot(Planeta[1], angtroy[0],1);
  
  	Planeta[2].Inicie(Xtroy+X0,Vtroy,m2,10);
	//-----------------------------------------
	vector3D Troy_noInerc;
	
	for (t=0; t<tmax; t+=dt){
		double theta = polar(Planeta[1]);
		Troy_noInerc = rot(Planeta[2], -theta,0);
	    	cout<<Troy_noInerc.x()<<"\t"<<Troy_noInerc.y()<<endl;
	    	
    		// Mover por PEFRL
    		for(i=0;i<N;i++) 
    			Planeta[i].Mueva_X(dt,Zeta); 
    			Newton.CalculeFuerzas(Planeta);
		for(i=0;i<N;i++) 
		   	Planeta[i].Mueva_V(dt,Coeficiente1);
    		for(i=0;i<N;i++)
    			Planeta[i].Mueva_X(dt,Chi);
    			Newton.CalculeFuerzas(Planeta);
    		for(i=0;i<N;i++) 
    			Planeta[i].Mueva_V(dt,Lambda);
    		for(i=0;i<N;i++)
    			Planeta[i].Mueva_X(dt,Coeficiente2);
    			Newton.CalculeFuerzas(Planeta);
    		for(i=0;i<N;i++)
    			Planeta[i].Mueva_V(dt,Lambda);
    		for(i=0;i<N;i++) 
    			Planeta[i].Mueva_X(dt,Chi);
    			Newton.CalculeFuerzas(Planeta);
    		for(i=0;i<N;i++)
    			Planeta[i].Mueva_V(dt,Coeficiente1);
    		for(i=0;i<N;i++) 
    			Planeta[i].Mueva_X(dt,Zeta); 
	}
	return 0;
}

