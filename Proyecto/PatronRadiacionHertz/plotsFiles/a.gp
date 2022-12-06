set term pdf
set out "Results/Comparison.pdf"
set xrange[55:100]
set g

J0=1
pi=3.14159265359

T=25.0
omega=2*pi/T
C=1/sqrt(2.0)
k=omega/C


alpha=0.5
s=1/sqrt(2*alpha)

p=J0/omega*s**3*(2*pi)**1.5

mu=2
eps=1
Z=sqrt(mu/eps)
A=Z*k*k*p/(4*pi)

phi = 0

t=70

abbs(x) = (x*x)**0.5
f(x) = A/x*abbs(cos(k*x-omega*t+phi)-sin(k*x-omega*t+phi)/(k*x))
g(x) = abbs(x-50)
h(x) = A/x*(cos(k*x-omega*t+phi)-sin(k*x-omega*t+phi)/(k*x))
plot "Datos/datazeL2.txt" u 1:(f(g($1))) w l title "teorico","Datos/datazeL2.txt" u 1:(abbs($2)) title "LB"
#plot "datazeL2.txt" u 1:((abbs($2))/(f(g($1)))) w l
#plot "Datos/datazeL2.txt" u 1:(h(g($1))) w l title "teorico","" u 1:2 title "LB"
