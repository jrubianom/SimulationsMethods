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

plotB:
	gnuplot plotsFiles/PolarB.gp
	xdg-open Results/PowerB.pdf

plotE:
	gnuplot plotsFiles/PolarE.gp
	xdg-open Results/PowerE.pdf

clean:
	rm *.txt *.x Datos/*.txt Results/*.pdf Datos/*.dat
