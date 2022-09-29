set g
set xl 'Temperatura k_B T'
set yl 'Presión'
set term pdf
set o 'datapunto_5i_presion.pdf'

plot 'data/punto_5i_presion.dat' w l notitle
