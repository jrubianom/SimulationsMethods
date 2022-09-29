set term pdf
set out "C.pdf"
set g
set logscale
f(x) = m*x+b
g(x) = m2*x+b2
fit f(x) "Tabla.txt" using (log($1)):(log($2)) via m,b
fit g(x) "Tabla.txt" using (log($1)):(log($3)) via m2,b2

ff(x) = exp(b)*x**m
gg(x) = exp(b2)*x**m2

Result = sprintf("T(k)=%.2f*k^{%.2f} \nt(k)=%.2f*k^{%.2f}",exp(b),m,exp(b2),m2)
set obj 1 rect from graph 0, 1 to graph 0.24, 0.86 fc rgb "white" front
set lab 1 Result at graph 0.01, 0.96 front

set xlabel "K"

plot "Tabla.txt" using 1:2 title "{/Symbol t medido} max" ,ff(x) title "{/Symbol t medido} Ajustado","Tabla.txt" using 1:3 title "t max" ,gg(x) title "t Ajustado"

set term pdf
set out "D.pdf"
unset obj 1
unset lab 1
unset logscale

tmin = system("awk 'NR == 1 {print $1}' Exp0.txt")
h(x,y) = (x-tmin)*(y**(-m2))
h2(x,y) = x*(y**(-m))

k1 = system("awk 'NR == 1 {print $1}' Tabla.txt")
k2 = system("awk 'NR == 2 {print $1}' Tabla.txt")
k3 = system("awk 'NR == 3 {print $1}' Tabla.txt")
k4 = system("awk 'NR == 4 {print $1}' Tabla.txt")
k5 = system("awk 'NR == 5 {print $1}' Tabla.txt")
k6 = system("awk 'NR == 6 {print $1}' Tabla.txt")

tmax = 3*system("awk 'NR == 1 {print $3}' Tabla.txt")+tmin
tmax = h(tmax,k1)
set xrange[0:tmax]
set xlabel sprintf(" t k^{%.2f}",-m2)
set ylabel sprintf("{/Symbol t} k^{%.2f}",-m)
plot "Exp0.txt" using (h($1,k1)):(h2($2,k1)) title "K1" w l,"Exp1.txt" using (h($1,k2)):(h2($2,k2)) title "K2" w l,"Exp2.txt" using (h($1,k3)):(h2($2,k3)) title "K3" w l,"Exp3.txt" using (h($1,k4)):(h2($2,k4)) title "K4" w l,"Exp4.txt" using (h($1,k5)):(h2($2,k5)) title "K5" w l,"Exp5.txt" using (h($1,k6)):(h2($2,k6)) title "K6" w l
