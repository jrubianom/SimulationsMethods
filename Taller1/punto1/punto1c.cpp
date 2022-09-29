#include <iostream>
#include <cmath>
#include <vector>

using namespace std;



//X(t)=(s(t),i(t),r(t))

double f(int M, double t, std::vector<double> X, double Gamma, double Beta){
	if (M==0){ 
		return -Beta*X[0]*X[1];
	}
	else if (M==1){
		return Beta*X[0]*X[1]-Gamma*X[1];
	}
	return 0;
}


void RK4dobleVector(double & t, std::vector<double> & Dx0, double Gamma, double Beta,  double dt){
	int N = Dx0.size(); std::vector<double> Dx1(N), Dx2(N), Dx3(N), Dx4(N), Aux(N);
	
	for(int i = 0; i < N; i++){
		Dx1[i] = dt*f(i,t, Dx0,Gamma,Beta);
	}
	for(int j = 0; j < N; j++){
		for (int jj = 0; jj < N; jj++){
			Aux[jj] = Dx0[jj]+Dx1[jj]/2;
		}
		Dx2[j] = dt*f(j,t+dt/2,Aux, Gamma, Beta);
	}
	for(int k = 0; k < N; k++){
		for (int kk = 0; kk < N; kk++){
			Aux[kk] = Dx0[kk]+Dx2[kk]/2;
		}
		Dx3[k] = dt*f(k,t+dt/2,Aux, Gamma, Beta);
	}
	for(int l = 0; l < N; l++){
		for (int ll = 0; ll< N; ll++){
			Aux[ll] = Dx0[ll]+Dx3[ll];
		}
		Dx4[l] = dt*f(l,t+dt,Aux, Gamma, Beta);
	}
	for (int n = 0; n < N; n++){
		Dx0[n] += (Dx1[n]+2*Dx2[n]+2*Dx3[n]+Dx4[n])/6;
	}
	t+= dt;
}

//obtiene Soo
double soo(double dt,double Gamma, double Beta){
	std::vector<double> X = {0.999, 0.001, 0}; 
	double t;
	double T = 1/Gamma;
	double Tmax = 1000*T;
	for(t=0; t<Tmax;){
		RK4dobleVector(t,X,Gamma, Beta, dt);
	}
	return X[0];
}

int main(){
	double dt = 0.01, dbeta = 0.001;
	double Gamma = 0.08;
	double s;
	for (double Beta = 0; Beta< 5*Gamma; Beta+= dbeta ){
		s = soo(dt,Gamma,Beta);
		cout<< Beta/Gamma <<"\t"<<s<<endl;
	}
	return 0;
}




