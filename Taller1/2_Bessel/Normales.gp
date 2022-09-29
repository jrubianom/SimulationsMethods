set term pdf
set g
set out "Modos0.pdf"
set xlabel "r"
set title "Solution R(r=1, {/Symbol l}=2.40483)"
set title font "Helvetica,14"
plot "Modo0.txt" using 1:3 w l title "Modo 0", "ModoExacto0.txt" using 1:2 w l title "Modo 0 Exacto"

set out "Modos1.pdf"
set title "Solution R(r=1, {/Symbol l}=5.52008)"
plot "Modo1.txt" using 1:3 w l title "Modo 1", "ModoExacto1.txt" using 1:2 w l title "Modo 1 Exacto"

set out "Modos2.pdf"
set title "Solution R(r=1, {/Symbol l}=8.65374)"
plot "Modo2.txt" using 1:3 w l title "Modo 2", "ModoExacto2.txt" using 1:2 w l title "Modo 2 Exacto"

set out "Modos3.pdf"
set title "Solution R(r=1, {/Symbol l}=11.7916)"
plot "Modo3.txt" using 1:3 w l title "Modo 3", "ModoExacto3.txt" using 1:2 w l title "Modo 3 Exacto"

set out "Modos4.pdf"
set title "Solution R(r=1, {/Symbol l} =14.9310)"
plot "Modo4.txt" using 1:3 w l title "Modo 4", "ModoExacto4.txt" using 1:2 w l title "Modo 4 Exacto"

set out "TodosModos.pdf"
set title "Modos normales experimentales R(r=1,{/Symbol l})"
set encoding iso_8859_1 
plot "Modo0.txt" using 1:3 w l title "{/Symbol l}_0 = 2.40483",\
"Modo1.txt" using 1:3 w l title "{/Symbol l}_1 = 5.52008",\
"Modo2.txt" using 1:3 w l title "{/Symbol l}_2 = 8.65374",\
"Modo3.txt" using 1:3 w l title "{/Symbol l}_3 = 11.7916",\
"Modo4.txt" using 1:3 w l title "{/Symbol l}_4 = 14.9310"