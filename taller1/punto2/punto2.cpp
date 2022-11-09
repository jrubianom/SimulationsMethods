#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double Tol = 1e-7;

//X(t)=(dR/dr(r),R(r))

double f(int M, double r, std::vector<double> X, double Lambda){
	if (M==0){ 
		return -X[1]*r*Lambda*Lambda;
	}
	else if (M==1){
		if (r<Tol){
			return 0;
		}
		else {
			return X[0]/r;
		}
	}
	return 0;
}


void RK4dobleVector(double & r, std::vector<double> & Dx0,  double dr, double Lambda){
	int N = Dx0.size(); std::vector<double> Dx1(N), Dx2(N), Dx3(N), Dx4(N), Aux(N);
	
	for(int i = 0; i < N; i++){
		Dx1[i] = dr*f(i,r, Dx0, Lambda);
	}
	for(int j = 0; j < N; j++){
		for (int jj = 0; jj < N; jj++){
			Aux[jj] = Dx0[jj]+Dx1[jj]/2;
		}
		Dx2[j] = dr*f(j,r+dr/2,Aux, Lambda);
	}
	for(int k = 0; k < N; k++){
		for (int kk = 0; kk < N; kk++){
			Aux[kk] = Dx0[kk]+Dx2[kk]/2;
		}
		Dx3[k] = dr*f(k,r+dr/2,Aux, Lambda);
	}
	for(int l = 0; l < N; l++){
		for (int ll = 0; ll< N; ll++){
			Aux[ll] = Dx0[ll]+Dx3[ll];
		}
		Dx4[l] = dr*f(l,r+dr,Aux, Lambda);
	}
	for (int n = 0; n < N; n++){
		Dx0[n] += (Dx1[n]+2*Dx2[n]+2*Dx3[n]+Dx4[n])/6;
	}
	r+= dr;
}




int main(){
	std::vector<double> X = {0, 1}; 
	double Lambda = 1;
	double r, dr = 0.01;
	double a = 1;
	for(r=0.01; r<10*a;){
		cout<<r<<"\t"<<X[1]<<endl;
		RK4dobleVector(r,X,dr,Lambda);
	}
	return 0;
}



