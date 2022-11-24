all:plot

%.x:%.cpp
	g++ $< -o $@
	./$@ $(p1) $(p2) $(N) $(t)

plot: main.x plot.gp
	gnuplot plot.gp
	xpdf Var.pdf
	rm main.x

clean:
	rm *.dat *.pdf *.log *.x
