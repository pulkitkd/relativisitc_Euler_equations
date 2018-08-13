set key font ",12"
set xtics font ",12"
set ytics font ",12"
set xlabel font ",12"
set ylabel font ",12"
set xran[0:1]

set xlabel 'x'
set ylabel 'Density'
p 'sol' u 1:2 t 'DG pol' w l lw 1, \
  'avg' u 1:2 t 'DG avg' w p pt 6, \
  'h' u 1:2 t 'DG avg' w p pt 2, \
  'sin-density-exact.dat' u 1:2 t 'Exact' w l lw 2

pause 1
reread
