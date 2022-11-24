all:plot

%.x:%.cpp
	g++ $< -o $@
	./$@

plot: main.x plot.gp
	gnuplot plot.gp
	xdg-open data.jpg

clean:
	rm *.dat *.pdf *.x *.log
