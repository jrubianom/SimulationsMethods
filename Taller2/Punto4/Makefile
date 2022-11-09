all: plot

%.x: Dielecttric_EMLBM.cpp
	g++ $< -o $@
	./$@ > data.dat

plot: Dielecttric_EMLBM.x plot.gp Dielecttric_EMLBM.cpp
	gnuplot plot.gp
	xdg-open data.jpg

clean:
	rm *.dat *.jpg *.x a.out
