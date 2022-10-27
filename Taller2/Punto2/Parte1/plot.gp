set term jpeg enhanced
set o "data.jpg"

set pm3d map
set size ratio 1
set xrange[0:200]
set yrange[0:200]
splot "data.dat"
