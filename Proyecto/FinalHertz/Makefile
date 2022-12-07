all:plot animation

Datos/dataBy.txt: main.cpp
	g++ $< -o naim.x -O3
	./naim.x


plot: plotsFiles/a.gp Datos/dataBy.txt
	gnuplot plotsFiles/FieldsInx.gp
	gnuplot plotsFiles/Polar.gp
	xdg-open Results/Comparison.pdf
	xdg-open Results/Power.pdf

replot:
	gnuplot plotsFiles/FieldsInx.gp
	gnuplot plotsFiles/Polar.gp
	xdg-open Results/Comparison.pdf
	xdg-open Results/Power.pdf

contour:
	gnuplot plotsFiles/BL.gp
	xpdf Results/Hertz.pdf

animation:
	gnuplot animationBy.gp
	gnuplot animationEz.gp
	gnuplot animationS.gp
	xdg-open FiledS.gif




clean:
	rm *.txt *.x Datos/*.txt Results/*.pdf
