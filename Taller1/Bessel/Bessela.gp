set term pdf
set out "Bessela.pdf"
set xlabel "r"
set title "Solution R"
set title font "Helvetica,14"
plot "Lambda1.txt" using 1:3 w l title "R"
