all:plot

naim.x: main.cpp
	g++ $< -o naim.x -O3
	./naim.x


plot: plotsFiles/a.gp naim.x
	gnuplot plotsFiles/FieldsInx.gp
	gnuplot plotsFiles/Polar.gp
	xdg-open Results/Comparison.pdf
	xdg-open Results/Power.pdf
	make contour

replot:
	gnuplot plotsFiles/FieldsInx.gp
	gnuplot plotsFiles/Polar.gp
	xdg-open Results/Comparison.pdf
	xdg-open Results/Power.pdf
	make contour

contour:
	gnuplot plotsFiles/BL.gp
	xpdf Results/Hertz.pdf

clean:
	rm *.txt *.x Datos/*.txt Results/*.pdf
