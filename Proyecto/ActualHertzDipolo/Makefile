all:plot

Datos/datazeL2.txt: main.cpp
	g++ $< -o naim.x -O3
	./naim.x


plot: plotsFiles/a.gp Datos/datazeL2.txt
	gnuplot plotsFiles/a.gp
	xdg-open Results/Comparison.pdf
	gnuplot plotFiles/c.gp
	xdg-open Results/Power.pdf

contour:
	gnuplot plotsFiles/BL.gp
	xpdf Results/Hertz.pdf

clean:
	rm *.txt *.x
