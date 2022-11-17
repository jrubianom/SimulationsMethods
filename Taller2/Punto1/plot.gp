set term pdf
set o "Var.pdf"
f(x)=m*x+b
fit f(x) "var.dat" u 1:2 via m,b
set xlabel "t"
set ylabel "Varianza"

Result = sprintf("{/Symbol s}^2=%.2f*t + {%.2f} ",m,b)
set obj 1 rect from graph 0, 1 to graph 0.24, 0.86 fc rgb "white" front
set lab 1 Result at graph 0.01, 0.96 front

plot "var.dat" u 1:2, f(x)
