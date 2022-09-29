set grid
set xl 'Tiempo'
set yl'Coordenada x'

set term pdf
set o 'data/punto_5a.pdf'

plot 'data/punto_5a.dat' w l notitle
