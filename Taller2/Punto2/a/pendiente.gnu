reset
clear
set term pdf
set output 'pendiente.pdf'
set xlabel "ix" font "TimesNewRoman,14"
set ylabel "iy" font "TimesNewRoman,14"
set grid
show grid


set dummy t
F(t) = a* t+ b
fit F(t) "Pendiente.dat" via a,b

plot "Pendiente.dat", F(t) w l
