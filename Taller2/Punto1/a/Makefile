all:plot

%.x:%.cpp
	g++ $< -o $@
	./$@ > var.dat

plot: main.x plot.gp
	gnuplot plot.gp
	xpdf Var.pdf

clean:
	rm *.dat *.pdf *.log *.x main.cpp
