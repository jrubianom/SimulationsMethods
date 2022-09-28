set term pdf
set out "Modosa.pdf"
set xlabel "r"
set title "Solution R"
set title font "Helvetica,14"
plot "Modosa.dat" using 1:2 w l title "R"