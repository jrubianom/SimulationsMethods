all:plot

%.x: %.cpp
	g++ $< -o $@
	./$@ > Datos.txt

plot: CunaNewton.x
	gnuplot gplot.gp
	xpdf datos.pdf
	gnuplot g2plot.gp
	xpdf Log.pdf
clean:
	rm ./a.out *.x
