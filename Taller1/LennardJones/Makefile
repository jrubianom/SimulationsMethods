all:

%.x: %.cpp Particulas.h
	g++ $< -o $@
	./$@

plot1: 1LJ.x 1LJ.cpp
	gnuplot Animacion1.txt
	gnuplot plot1.gp
	xpdf Energia.pdf
	xpdf xvst.pdf
	xdg-open UnaParticula.gif

plot2:LJ.x  Particulas.h
	gnuplot 2DGas.txt
	xdg-open 2DGas.gif
	gnuplot plot2.gp
	xpdf Energia2.pdf
	xpdf yprom.pdf

cleandata:
	rm *.txt *.pdf *.x a.out
clean:
	rm *.x a.out
