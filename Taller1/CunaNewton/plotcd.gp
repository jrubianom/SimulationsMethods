set term pdf
set out "C-D.pdf"
f(x) = m*x+b
g(x) = m2*x+b2
fit f(x) "Tabla.txt" using (log($1)):(log($2)) via m,b
fit g(x) "Tabla.txt" using (log($1)):(log($3)) via m2,b2
plot "Tabla.txt" using (log($1)):(log($2)) title "{/Symbol t} max" ,f(x) title "fit {/Symbol t}","Tabla.txt" using (log($1)):(log($3)) title "t max" ,g(x) title "fit t"

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
set xlabel "t"
set ylabel "{/Symbol t}"
plot "Exp0.txt" using (h($1,k1)):(h2($2,k1)) title "K1" w l,"Exp1.txt" using (h($1,k2)):(h2($2,k2)) title "K2" w l,"Exp2.txt" using (h($1,k3)):(h2($2,k3)) title "K3" w l,"Exp3.txt" using (h($1,k4)):(h2($2,k4)) title "K4" w l,"Exp4.txt" using (h($1,k5)):(h2($2,k5)) title "K5" w l,"Exp5.txt" using (h($1,k6)):(h2($2,k6)) title "K6" w l
