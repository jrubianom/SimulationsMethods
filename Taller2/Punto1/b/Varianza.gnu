reset
clear
set terminal jpeg enhanced
set output 'VarianzaVsT.jpg'
set xlabel "t" font "TimesNewRoman,14"
set ylabel "Varianza" font "TimesNewRoman,14"
set grid
show grid


set dummy t
F(t) = a* t+ b
fit F(t) "VarianzaVsT.txt" via a,b

plot "VarianzaVsT.txt", F(t) w l
