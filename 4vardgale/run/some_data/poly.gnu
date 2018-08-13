set key font ",12"
set xtics font ",12"
set ytics font ",12"
set xlabel font ",12"
set ylabel font ",12"

set xlabel 'x'
set ylabel 'Density'
p 'sol' u 1:2 t 'After Transfer' w l lw 1, \
  'solo' u 1:2 t 'Before transfer' w l lw 1, \
  'avg' u 1:2 t 'DG Avg' w p pt 6, \
  'poly.avg' u 1:2 t 'Exact' w p pt 2, \
  'poly.dat' u 1:2 t 'Exact' w l lw 0.5

pause 1
reread
