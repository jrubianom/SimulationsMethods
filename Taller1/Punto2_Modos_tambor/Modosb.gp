set term pdf
set out "Modosb.pdf"
set xlabel "{/Symbol l}"
set title "Eigenvalues of R"
set title font "Helvetica,14"
plot "Modosb.dat" using 1:2 w l title "{/Symbol l} axis","Modosb.dat" using 1:3 w l title "R(1,{/Symbol l})"