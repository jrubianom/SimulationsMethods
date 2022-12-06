set term pdf
set out "Results/Power.pdf"
set g

set mapping spherical
splot "Datos/PowerDistri.txt" using 1:2:3 w l
