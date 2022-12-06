set term pdf
set out "Results/Comparison.pdf"
set xrange[55:100]
set g

plot "Datos/dataBy.txt" u 1:3 w l title "teorico","" u 1:2 title "LB"
set yrange[-0.1:0.1]
plot "Datos/dataEz.txt" u 1:3 w l title "teorico","" u 1:2 title "LB"
