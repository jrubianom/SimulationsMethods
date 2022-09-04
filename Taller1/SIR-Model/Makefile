all: plot1 plot2

%.x: %.cpp
	g++ $< -o $@
%.txt: %.x
	./$< > $@

%.pdf: %.txt %.gp
	cp $< datos.txt
	gnuplot $(filter %.gp,$^)

plot1: SIRa.pdf
	xpdf $<

plot2: SIRc.pdf
	xpdf $<

clean:
	rm *txt *x datos*
