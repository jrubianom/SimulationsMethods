reset
clear
set terminal jpeg enhanced
set output 'perfil100.jpg'
set xlabel "ix" font "TimesNewRoman,14"
set ylabel "rho" font "TimesNewRoman,14"
set grid
show grid


plot "Perfil100.dat" w l

reset
clear
set terminal jpeg enhanced
set output 'penrfil0.jpg'
set xlabel "ix" font "TimesNewRoman,14"
set ylabel "rho" font "TimesNewRoman,14"
set grid
show grid


plot "Perfil0.dat"  w l
