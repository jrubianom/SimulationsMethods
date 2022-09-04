set term pdf
set out "Log.pdf"
set logscale x
set logscale y
plot "Datos.txt" using 1:2 title "{/Symbol L} max","Datos.txt" using 1:3 title "t max"
