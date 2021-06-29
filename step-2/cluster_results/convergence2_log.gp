set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 2 pt 7 ps 0.5
set style line 3 lc rgb "green" lt 1 lw 2 pt 7 ps 0.5
set style line 4 lc rgb "orange" lt 1 lw 2 pt 7 ps 0.5

set title 'Refinaments Vs. Total Mean Error'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set terminal jpeg enhanced
set output 'refs_error_log.jpg'
set logscale x
set logscale y

plot 'refs.txt' u 1:2 w l ls 1 t 'solver 8', 'refs.txt' u 1:3 w l ls 2 t 'solver 9', 'refs.txt' u 1:4 w l ls 3 t 'solver 11', 'refs.txt' u 1:5 w l ls 4 t 'solver 12' 

unset title
unset ylabel

set title 'Refinaments Vs. Execution Time'
set ylabel 'Time [min]'
set output 'refs_time_log.jpg'

plot 'refs.txt' u 1:6 w l ls 1 t 'solver 8', 'refs.txt' u 1:7 w l ls 2 t 'solver 9', 'refs.txt' u 1:8 w l ls 3 t 'solver 11', 'refs.txt' u 1:9 w l ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Refinaments Vs. Iterations'
set ylabel 'Iterations'
set output 'refs_iter_log.jpg'

plot 'refs.txt' u 1:10 w l ls 1 t 'solver 8', 'refs.txt' u 1:11 w l ls 2 t 'solver 9', 'refs.txt' u 1:12 w l ls 3 t 'solver 11', 'refs.txt' u 1:13 w l ls 4 t 'solver 12'

unset title
unset ylabel
unset xlabel

set title 'Order Vs. Total Mean Error'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set output 'order_error_log.jpg'

plot 'order.txt' u 1:2 w l ls 1 t 'solver 8', 'order.txt' u 1:3 w l ls 2 t 'solver 9', 'order.txt' u 1:4 w l ls 3 t 'solver 11', 'order.txt' u 1:5 w l ls 4 t 'solver 12' 

unset title
unset ylabel

set title 'Order Vs. Execution Time'
set ylabel 'Time [min]'
set output 'order_time_log.jpg'

plot 'order.txt' u 1:6 w l ls 1 t 'solver 8', 'order.txt' u 1:7 w l ls 2 t 'solver 9', 'order.txt' u 1:8 w l ls 3 t 'solver 11', 'order.txt' u 1:9 w l ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Order Vs. Iterations'
set ylabel 'Iterations'
set output 'order_iter_log.jpg'

plot 'order.txt' u 1:10 w l ls 1 t 'solver 8', 'order.txt' u 1:11 w l ls 2 t 'solver 9', 'order.txt' u 1:12 w l ls 3 t 'solver 11', 'order.txt' u 1:13 w l ls 4 t 'solver 12'

