reset
clear
set term png
set output 'SIR.png'
set xlabel "Tiempo" font "TimesNewRoman,14"
set ylabel "Fracción de la población" font "TimesNewRoman,14"
set grid
show grid

plot 'SIR.dat' u 1:2 w l lw 2 title "S(t)", 'SIR.dat' u 1:3 w l lw  2 title "I (t)", 'SIR.dat' u 1:4 w l lw 2 title "R(t)"
