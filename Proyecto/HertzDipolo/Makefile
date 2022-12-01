all:plot

datazeL2.txt: main.cpp
	g++ $< -o naim.x
	./naim.x


plot: a.gp datazeL2.txt
	gnuplot a.gp
	xdg-open Comparison.pdf

contour:
	gnuplot BL.gp
	xpdf Hertz.pdf
