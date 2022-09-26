set term pdf
set out "A.pdf"
set g
set xlabel "t [s]"
tmin = system("awk 'NR == 1 {print $1}' Exp0.txt")
tmax = 3*system("awk 'NR == 1 {print $3}' Tabla.txt")+tmin
set xrange[tmin:tmax]
set ylabel "{/Symbol t} [Dina cm]"
plot "Exp0.txt" using 1:2 w l title "K1","Exp1.txt" using 1:2 w l title "K2","Exp2.txt" using 1:2 w l title "K3","Exp3.txt" using 1:2 w l title "K4","Exp4.txt" using 1:2 w l title "K5","Exp5.txt" using 1:2 w l title "K6"
