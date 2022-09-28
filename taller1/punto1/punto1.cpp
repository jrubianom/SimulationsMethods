#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


const double Beta  = 0.35;
const double Gamma = 0.08;
const double R0 = Beta/Gamma;

//X(t)=(s(t),i(t),r(t))

double f(int M, double t, std::vector<double> X){
	if (M==0){ 
		return -Beta*X[0]*X[1];
	}
	else if (M==1){
		return Beta*X[0]*X[1]-Gamma*X[1];
	}
	return 0;
}
// Función que selecciona las ecuaciones diferenciales a resolver M=0 para s(t) y M=1 para i(t); r(t)=1-s(t)-i(t).

// Implementación de RungeKutta

void RK4dobleVector(double & t, std::vector<double> & Dx0,  double dt){
	int N = Dx0.size(); std::vector<double> Dx1(N), Dx2(N), Dx3(N), Dx4(N), Aux(N);
	
	for(int i = 0; i < N; i++){
		Dx1[i] = dt*f(i,t, Dx0);
	}
	for(int j = 0; j < N; j++){
		for (int jj = 0; jj < N; jj++){
			Aux[jj] = Dx0[jj]+Dx1[jj]/2;
		}
		Dx2[j] = dt*f(j,t+dt/2,Aux);
	}
	for(int k = 0; k < N; k++){
		for (int kk = 0; kk < N; kk++){
			Aux[kk] = Dx0[kk]+Dx2[kk]/2;
		}
		Dx3[k] = dt*f(k,t+dt/2,Aux);
	}
	for(int l = 0; l < N; l++){
		for (int ll = 0; ll< N; ll++){
			Aux[ll] = Dx0[ll]+Dx3[ll];
		}
		Dx4[l] = dt*f(l,t+dt,Aux);
	}
	for (int n = 0; n < N; n++){
		Dx0[n] += (Dx1[n]+2*Dx2[n]+2*Dx3[n]+Dx4[n])/6;
	}
	t+= dt;
}




int main(){
	//vector con los valores iniciales
	//----------------------(s(0) , i(0) ,r())
	std::vector<double> X = {0.999, 0.001, 0}; 
	double t; double dt = 0.01;
	double T = 1/Gamma; //tiempo característico de la recuperación.
	double Tmax = 15*T; //tiempo máximo de la sumilación
	for(t=0; t<Tmax;){
		// ds/dt = -beta*s*i
		cout<<t<<"\t"<<X[0]<<"\t"<<X[1]<<"\t"<<1-X[0]-X[1]<<endl;
		RK4dobleVector(t,X,dt);
	}
	return 0;
}




