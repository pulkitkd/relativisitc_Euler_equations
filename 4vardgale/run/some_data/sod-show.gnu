set key font ",12"
set xtics font ",12"
set ytics font ",12"
set xlabel font ",12"
set ylabel font ",12"

set xlabel 'x'
set ylabel 'Density'
p 'sod.sol' u 1:2 t 'Refinment' w l lw 2, \
  'sod-no-ref.dat' u 1:2 t 'No Refinment' w l lw 2, \
  'sod-no-ref.avg' u 1:2 t 'No Ref Average' w p pt 6, \
  'sod.avg' u 1:2 t 'Ref Avg' w p pt 6, \
  'sod.dat' u 1:2 t 'Exact' w l lw 2
