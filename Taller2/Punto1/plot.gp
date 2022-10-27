set term pdf
set o "Var.pdf"
f(x)=m*x+b
fit f(x) "var.dat" u 1:2 via m,b
plot "var.dat" u 1:2, f(x)
