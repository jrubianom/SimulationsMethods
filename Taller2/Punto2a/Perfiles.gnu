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
set output 'penrfil60.jpg'
set xlabel "ix" font "TimesNewRoman,14"
set ylabel "rho" font "TimesNewRoman,14"
set grid
show grid


plot "Perfil60.dat"  w l
