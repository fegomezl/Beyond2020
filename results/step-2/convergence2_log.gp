set g
set term pdf
set key opaque
set key t r

set ls 1 lc rgb 'blue' lt 1 pt 1 lw 1 ps 1
set ls 2 lc rgb 'red' lt 2 pt 2 lw 1 ps 1
set ls 3 lc rgb 'black' lt 6 pt 6 lw 1 ps 1
set ls 4 lc rgb 'green' lt 6 pt 4 lw 1 ps 1

set title 'Refinaments Vs. Total Mean Error'
set xlabel 'Mesh size (u.a)'
set ylabel 'Total Mean Error'

set output 'refs_error_log.pdf'
set logscale x
set logscale y

plot 'refs.txt' u 1:2  ls 1 t 'solver 8', 'refs.txt' u 1:3  ls 2 t 'solver 9', 'refs.txt' u 1:4  ls 3 t 'solver 11', 'refs.txt' u 1:5  ls 4 t 'solver 12' 

unset title
unset ylabel

set title 'Refinaments Vs. Execution Time'
set ylabel 'Time [min]'
set output 'refs_time_log.pdf'

plot 'refs.txt' u 1:6  ls 1 t 'solver 8', 'refs.txt' u 1:7  ls 2 t 'solver 9', 'refs.txt' u 1:8  ls 3 t 'solver 11', 'refs.txt' u 1:9  ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Refinaments Vs. Iterations'
set ylabel 'Iterations'
set output 'refs_iter_log.pdf'

plot 'refs.txt' u 1:10  ls 1 t 'solver 8', 'refs.txt' u 1:11  ls 2 t 'solver 9', 'refs.txt' u 1:12  ls 3 t 'solver 11', 'refs.txt' u 1:13  ls 4 t 'solver 12'

unset title
unset ylabel
unset xlabel
unset logscale x
set xtics 1

set title 'Polinomyal order Vs. Total Mean Error'
set xlabel 'Order'
set ylabel 'Total Mean Error'
set output 'order_error_log.pdf'

plot 'order.txt' u 1:2  ls 1 t 'solver 8', 'order.txt' u 1:3 ls 2 t 'solver 9', 'order.txt' u 1:4  ls 3 t 'solver 11', 'order.txt' u 1:5  ls 4 t 'solver 12' 

unset title
unset ylabel

set title 'Order Vs. Execution Time'
set ylabel 'Time [min]'
set output 'order_time_log.pdf'

plot 'order.txt' u 1:6  ls 1 t 'solver 8', 'order.txt' u 1:7  ls 2 t 'solver 9', 'order.txt' u 1:8  ls 3 t 'solver 11', 'order.txt' u 1:9  ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Order Vs. Iterations'
set ylabel 'Iterations'
set output 'order_iter_log.pdf'

plot 'order.txt' u 1:10  ls 1 t 'solver 8', 'order.txt' u 1:11  ls 2 t 'solver 9', 'order.txt' u 1:12  ls 3 t 'solver 11', 'order.txt' u 1:13  ls 4 t 'solver 12'

