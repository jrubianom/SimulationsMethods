set term pdf
set out "Modos0.pdf"
set xlabel "r"
set title "Solution R"
set title font "Helvetica,14"
plot "Modo0.txt" using 1:3 w l title "Modo 0", "ModoExacto0.txt" using 1:2 w l title "Modo 0 Exacto"

set out "Modos1.pdf"
plot "Modo1.txt" using 1:3 w l title "Modo 1", "ModoExacto1.txt" using 1:2 w l title "Modo 1 Exacto"

set out "Modos2.pdf"
plot "Modo2.txt" using 1:3 w l title "Modo 2", "ModoExacto2.txt" using 1:2 w l title "Modo 2 Exacto"

set out "Modos3.pdf"
plot "Modo3.txt" using 1:3 w l title "Modo 3", "ModoExacto3.txt" using 1:2 w l title "Modo 3 Exacto"

set out "Modos4.pdf"
plot "Modo4.txt" using 1:3 w l title "Modo 4", "ModoExacto4.txt" using 1:2 w l title "Modo 4 Exacto"
