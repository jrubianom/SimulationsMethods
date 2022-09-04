all: plot

%.x: %.cpp
	g++ $< -o $@
%.txt: %.x
	./$< > $@

%.pdf: %.txt %.gp
	cp $< datos.txt
	gnuplot $(filter %.gp,$^)

plot1: Bessela.pdf
	xpdf $<

plot2: Besselb.pdf
	xpdf $<

table: EigenBessel.cpp
	g++ $<
	./a.out
clean:
	rm *txt *x *pdf a.out
