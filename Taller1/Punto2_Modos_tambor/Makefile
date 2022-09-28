all: plot1 plot2 table

%.x: %.cpp
	g++ $< -o $@

%.dat: %.x
	./$< > $@

%.pdf: %.dat %.gp
	gnuplot $(filter %.gp,$^)

plot1: Modosa.pdf
	xpdf $<

plot2: Modosb.pdf
	xpdf $<

table: EigenModos.cpp
	g++ $<
	./a.out
	gnuplot Normales.gp
	xpdf Modos0.pdf
	xpdf Modos1.pdf
	xpdf Modos2.pdf
	xpdf Modos3.pdf
	xpdf Modos4.pdf

clean:
	rm *.dat *.x *.pdf a.out