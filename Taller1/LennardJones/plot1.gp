set term pdf
set out "Energia.pdf"
set xlabel "t"
set ylabel "E m"
set yrange[-510:-490]
plot "datos1.txt" u 1:3 w l title "Energia total"

set out "xvst.pdf"
unset yrange
set xlabel "t"
set ylabel "x"
set title "Osiclación al rededor \ndel punto de equilibrio"
plot "datos1.txt" u 1:2 w l
