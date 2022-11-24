all:plot

%.x:%.cpp
	g++ $< -o $@
	./$@

plot: main.x animation.gp
	gnuplot animation.gp
	xdg-open Mirror.gif

clean:
	rm *.dat *.pdf *.x *.log
