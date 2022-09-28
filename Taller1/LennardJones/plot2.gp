set term pdf
set out "Energia2.pdf"
set g
set logscale
set xlabel "t"
set ylabel "E m"
plot "data.txt" u 1:2 w l title "Energia total"

set out "yprom.pdf"
unset logscale
set xlabel "t"
set ylabel "yprom"
set title "Altura Promedio"
plot "data.txt" u 1:3 w l
