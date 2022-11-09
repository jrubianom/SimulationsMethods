set term jpeg enhanced
set o "data.jpg"

set size ratio 1
plot "data.dat" u 1:2 w l
