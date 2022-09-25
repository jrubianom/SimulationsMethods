set term pdf
set out "C-D.pdf"
f(x) = m*x+b
g(x) = m2*x+b2
fit f(x) "Tabla.txt" using (log($1)):(log($2)) via m,b
fit g(x) "Tabla.txt" using (log($1)):(log($3)) via m2,b2
plot "Tabla.txt" using (log($1)):(log($2)) title "{/Symbol t} max" ,f(x) title "fit {/Symbol t}","Tabla.txt" using (log($1)):(log($3)) title "t max" ,g(x) title "fit t"

h(x,y) = (x-0.17456)*(y**(-m2))
h2(x,y) = x*(y**(-m))
h3(x,y) = x*(y**(-m2))

k1 = system("awk 'NR == 1 {print $1}' Tabla.txt")
k2 = system("awk 'NR == 2 {print $1}' Tabla.txt")
k3 = system("awk 'NR == 3 {print $1}' Tabla.txt")
k4 = system("awk 'NR == 4 {print $1}' Tabla.txt")
k5 = system("awk 'NR == 5 {print $1}' Tabla.txt")
k6 = system("awk 'NR == 6 {print $1}' Tabla.txt")

TT1= system("awk 'NR == 1 {print $2}' Tabla.txt")
TT2= system("awk 'NR == 2 {print $2}' Tabla.txt")
TT3= system("awk 'NR == 3 {print $2}' Tabla.txt")
TT4= system("awk 'NR == 4 {print $2}' Tabla.txt")
TT5= system("awk 'NR == 5 {print $2}' Tabla.txt")
TT6= system("awk 'NR == 6 {print $2}' Tabla.txt")

tt1= system("awk 'NR == 1 {print $3}' Tabla.txt")
tt2= system("awk 'NR == 2 {print $3}' Tabla.txt")
tt3= system("awk 'NR == 3 {print $3}' Tabla.txt")
tt4= system("awk 'NR == 4 {print $3}' Tabla.txt")
tt5= system("awk 'NR == 5 {print $3}' Tabla.txt")
tt6= system("awk 'NR == 6 {print $3}' Tabla.txt")

T1=h2(TT1,k1)
T2=h2(TT2,k2)
T3=h2(TT3,k3)
T4=h2(TT4,k4)
T5=h2(TT5,k5)
T6=h2(TT6,k6)
t1=h3(tt1,k1)
t2=h3(tt2,k2)
t3=h3(tt3,k3)
t4=h3(tt4,k4)
t5=h3(tt5,k5)
t6=h3(tt6,k6)

set logscale x
set xlabel "t"
set ylabel "{/Symbol t}"
plot "Exp0.txt" using (h($1,k1)/t1):(h2($2,k1)/T1) title "K1" w l,"Exp1.txt" using (h($1,k2)/t2):(h2($2,k2)/T2) title "K2" w l,"Exp2.txt" using (h($1,k3)/t3):(h2($2,k3)/T3) title "K3" w l,"Exp3.txt" using (h($1,k4)/t4):(h2($2,k4)/T4) title "K4" w l,"Exp4.txt" using (h($1,k5)/t5):(h2($2,k5)/T5) title "K5" w l,"Exp5.txt" using (h($1,k6)/t6):(h2($2,k6)/T6) title "K6" w l
