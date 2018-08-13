set xtics font ",20"
set ytics font ",20"
set xlabel font ",24"
set ylabel font ",24"
set xran[0:1]

set xlabel 'x'
set ylabel 'Density'
p 'sol' u 1:2 t 'DG pol' w l lw 2, \
  'avg' u 1:2 t 'DG avg' w p pt 6, \
  'sod.dat' u 1:2 t 'Exact' w l lw 2
