all: plot1 plot2 table

%.x: %.cpp
	g++ $< -o $@

%.txt: %.x
	./$< > $@

%.pdf: %.txt %.gp
	gnuplot $(filter %.gp,$^)

plot1: Bessela.pdf
#	xpdf $< &

plot2: Besselb.pdf
#	xpdf $< &

#Puntos c y d
table: EigenBessel.cpp
	g++ $<
	./a.out > Tabla.dat
	gnuplot Normales.gp
#	xpdf Modos0.pdf
#	xpdf Modos1.pdf
#	xpdf Modos2.pdf
#	xpdf Modos3.pdf
#	xpdf Modos4.pdf
#	xpdf TodosModos.pdf

clean:
	rm *.txt *.x a.out

cleanall:
	rm *.txt *.x a.out *.pdf *.dat
