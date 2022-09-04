set term pdf
set out "Bessela.pdf"
set xlabel "r"
set title "Solution R"
set title font "Helvetica,14"
plot "datos.txt" using 1:3 w l title "R"
