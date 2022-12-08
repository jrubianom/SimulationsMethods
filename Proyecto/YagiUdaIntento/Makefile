all:plot

%.x: main.cpp
	g++ $< -o main.x -O3
	./main.x


plot: plotsFiles/a.gp main.x
	gnuplot plotsFiles/Polar.gp
	xdg-open Results/Power.pdf

replot:
	gnuplot plotsFiles/Polar.gp
	xdg-open Results/Power.pdf

clean:
	rm *.txt *.x Datos/*.txt Results/*.pdf Datos/*.dat
