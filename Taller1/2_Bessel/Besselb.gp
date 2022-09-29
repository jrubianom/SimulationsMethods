set term pdf
set out "Besselb.pdf"
set xlabel "{/Symbol l}"
set ylabel "R"
set title "Eigenvalues of R"
set title font "Helvetica,14"
set mxtics 4
set xtics 2
set mytics 1
set ytics 0.2
set grid mxtics xtics
set grid mytics ytics
plot "Besselb.txt" using 1:2 w l title "{/Symbol l} axis","Besselb.txt" using 1:3 w l title "R(1,{/Symbol l})"
