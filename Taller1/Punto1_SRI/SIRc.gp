set term pdf
set out "SIRc.pdf"
set xlabel "beta/gamma"
set title "SIR model"
set title font "Helvetica,14"
plot "datos.txt" using 1:2 w l title "s"
