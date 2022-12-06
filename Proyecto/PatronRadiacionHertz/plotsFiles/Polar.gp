set term pdf
set out "Results/Power.pdf"
set polar
set g

f(x)=x*0.35/0.4

plot "Datos/EPlane.txt" u 1:2 title "S_E LB","Datos/TeoEPlane.txt" u 1:(f($2)) w l title "S_E Teorico","Datos/BPlane.txt" u 1:2 title "S_B LB","Datos/TeoBPlane.txt" u 1:(f($2)) w l title "S_B Teorico"
