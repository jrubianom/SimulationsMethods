set term pdf
set out "Results/Power.pdf"
set polar

set g

f(x)=x*17.5/20

plot "Datos/BPlane.txt" u 1:2 w l title "S_B LB"
plot "Datos/EPlane.txt" u 1:2 w l title "S_E LB"
#plot "Datos/BPlane.txt" u 1:2 title "S_B LB"
