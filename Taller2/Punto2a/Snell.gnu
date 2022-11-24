reset
clear
set pm3d map
set size ratio 1

set terminal jpeg enhanced
set output "Ondas.jpg"
splot "Ondas.dat"
