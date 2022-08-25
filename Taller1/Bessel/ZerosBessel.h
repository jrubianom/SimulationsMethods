#include <iostream>
#include <cmath>

double ErrMax = 1e-8;

double f(double alpha, double  x, double t){
  return std::cos(alpha*t-x*std::sin(t));
}

double IntegralPorSimpson(double alpha, double x,
                           double a, double b, int n){
  double Integral = 0.0;
  double h;
  n *= 2; h = (b-a)/n;
  Integral = f(alpha,x,a) + f(alpha,x,b) + 4*f(alpha,x,a + h);
  for(int i = 1; i < n/2; i++){
    Integral += 2*f(alpha,x,a + 2*i*h) + 4*f(alpha,x,a + (2*i+1)*h);
  }
  return h/3*Integral;
}

double Bessel_alpha(double alpha, double x,int n){
  double a = 0, b = M_PI;
  return IntegralPorSimpson(alpha,x,a,b,n)*1.0/M_PI;
}

double CerosPorBiseccion(double alpha, double a, double b,int n){
  double  m, fa, fm;
  fa=Bessel_alpha(alpha,a,n);
  while(b-a >= ErrMax){
    m = (b+a)/2; fm=Bessel_alpha(alpha,m,n);
    if(fa*fm>0)
    {a=m; fa=fm;}
    else
      b=m;
  }
  return (a+b)/2;
}

double ZeroBessel(double alpha,double a,double b,int n = 50){
  return CerosPorBiseccion(alpha,a,b,n);
}
