set terminal x11 noenhanced
set xlabel 'x'
set ylabel 'Density'
set xran[0:1]
files = system("ls -1 -v *.sol")
labels = system("ls -1 -v *.sol | sed -e 's/.sol//'")
do for [i=1:words(files)-1] { plot sprintf('%d.sol', i) u 1:2 title word(labels,i) w l lw 1, sprintf('%d.avg', i) u 1:2 title 'AVG' w p pt 6; pause 0.01 }
