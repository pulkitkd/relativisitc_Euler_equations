set key font ",12"
set xtics font ",12"
set ytics font ",12"
set xlabel font ",12"
set ylabel font ",12"

set xlabel 'x'
set ylabel 'Density'
p 'lowd.sol' u 1:2 t 'DG pol' w l lw 1, \
  'lowd.avg' u 1:2 t 'DG avg' w p pt 6, \
  'lowd-no-ref.sol' u 1:2 t 'DG pol' w l lw 1, \
  'lowd-no-ref.avg' u 1:2 t 'DG avg' w p pt 6, \
  'lowd-exact.avg' u 1:2 t 'DG avg' w p pt 2, \
  'lowd.dat' u 1:2 t 'Exact' w l lw 0.5
