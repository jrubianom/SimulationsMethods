set term pdf
set out "Bessela.pdf"
set xlabel "r"
set ylabel "R"
set title "Solution R(r)"
set title font "Helvetica,14"
set g
plot "Lambda1.txt" using 1:3 w l title "R(r)"
