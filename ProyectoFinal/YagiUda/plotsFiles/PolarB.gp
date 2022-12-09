set term pdf
set out "Results/PowerB.pdf"
set polar

set g

plot "Datos/BPlane.txt" u 1:2 title "S_B LB"
