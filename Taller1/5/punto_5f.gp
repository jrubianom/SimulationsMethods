set g
set xl 'Tiempo'
set yl 'Coordenada Y'

set term pdf
set o 'data/punto_5f.pdf'

plot for [i=2:26] 'data/punto_5f.dat' u 1:i w l notitle
