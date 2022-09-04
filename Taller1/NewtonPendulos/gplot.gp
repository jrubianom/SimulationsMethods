set term pdf
set out "datos.pdf"
set xrange[0.174:0.178]
plot "Exp0.txt" using 1:2 w l title "K1","Exp1.txt" using 1:2 w l title "K2","Exp2.txt" using 1:2 w l title "K3","Exp3.txt" using 1:2 w l title "K4","Exp4.txt" using 1:2 w l title "K5","Exp5.txt" using 1:2 w l title "K6"
