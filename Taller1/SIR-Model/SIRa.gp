set term pdf
set out "SIRa.pdf"
set xlabel "time"
set title "SIR model"
set title font "Helvetica,14"
plot "datos.txt" using 1:2 w l title "s", "datos.txt" using 1:3 w l title "i","datos.txt" using 1:4 w l title "r"
