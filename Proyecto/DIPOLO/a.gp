set term pdf
set out "Comparison.pdf"
set xrange[58:100]
set g

J0=1
pi=3.14159265359

T=25.0
omega=2*pi/T
C=1/sqrt(2.0)
k=omega/C


alpha=0.25
s=1/sqrt(2*alpha)

p=J0/omega*s**3*(2*pi)**1.5

mu=2
eps=1
Z=sqrt(mu/eps)
A=Z*k*k*p/(4*pi)



t=70

abbs(x) = (x*x)**0.5
f(x) = A/x*abbs(cos(k*x-omega*t)-sin(k*x-omega*t)/(k*x))
g(x) = abbs(x-50)
plot "datazeL2.txt" u 1:(f(g($1))) w l title "teorico","datazeL2.txt" u 1:(abbs($2)) title "LB"
#plot "datazeL2.txt" u 1:((abbs($2))/(f(g($1)))) w l
