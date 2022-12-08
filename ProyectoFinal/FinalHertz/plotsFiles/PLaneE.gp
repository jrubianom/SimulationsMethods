set term pdf
set out "Results/PowerPlaneE.pdf"
set polar
set g

plot "Datos/EPlane.txt" u 1:(2*$2/25.0) title "S_E LB","Datos/TeoEPlane.txt" u 1:2 w l title "S_E Teorico"
