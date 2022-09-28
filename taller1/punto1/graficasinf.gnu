reset
clear
set term png
set output 's_inf.png'
set xlabel "{/Symbol b}/{/Symbol g}" font "TimesNewRoman,14"
set ylabel "s({/Symbol \245})" font "TimesNewRoman,14"
set grid
show grid

plot 'sinf.dat' w l lw 2 title "s({/Symbol \245})"
