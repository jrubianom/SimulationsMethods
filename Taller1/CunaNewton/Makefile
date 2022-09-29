all:plota plotcd

%.x: %.cpp
	g++ $< -o $@
	./$@ > Tabla.txt

plota: CunaNewton.x
	gnuplot plota.gp
	xpdf A.pdf

plotcd:
	gnuplot plotcd.gp
	xpdf C.pdf
	xpdf D.pdf

Animacion: CunaNewton.x
	gnuplot AnimacionData.txt

clean:
	rm ./a.out *.x *.txt *.pdf *.gif *.log
