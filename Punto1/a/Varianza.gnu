reset
clear
set term pdf
set output 'VarianzaVsT.pdf'
set xlabel "t" font "TimesNewRoman,14"
set ylabel "Varianza" font "TimesNewRoman,14"
set grid
show grid


set dummy t
F(t) = a* t+ b
fit F(t) "VarianzaVsT.txt" via a,b

plot "VarianzaVsT.txt", F(t) w l
