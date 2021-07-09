set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 pt 1
set style line 2 lc rgb "red" lt 1 pt 2
set style line 3 lc rgb "green" lt 1 pt 3
set style line 4 lc rgb "black" lt 1 pt 6

set title 'Refinaments Vs. Total Mean Error'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set term pdf font "arial"
set output 'refs_error.pdf'

plot 'refs.txt' u 1:2 ls 1 t 'solver 8', 'refs.txt' u 1:3 ls 2 t 'solver 9', 'refs.txt' u 1:4 ls 3 t 'solver 11', 'refs.txt' u 1:5 ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Refinaments Vs. Execution Time'
set ylabel 'Time [min]'
set output 'refs_time.pdf'

plot 'refs.txt' u 1:6 ls 1 t 'solver 8', 'refs.txt' u 1:7 ls 2 t 'solver 9', 'refs.txt' u 1:8 ls 3 t 'solver 11', 'refs.txt' u 1:9 ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Refinaments Vs. Iterations'
set ylabel 'Iterations'
set output 'refs_iter.pdf'

plot 'refs.txt' u 1:10 ls 1 t 'solver 8', 'refs.txt' u 1:11 ls 2 t 'solver 9', 'refs.txt' u 1:12 ls 3 t 'solver 11', 'refs.txt' u 1:13 ls 4 t 'solver 12'

unset title
unset ylabel
unset xlabel

set title 'Order Vs. Total Mean Error'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set output 'order_error.pdf'

plot 'order.txt' u 1:2 ls 1 t 'solver 8', 'order.txt' u 1:3 ls 2 t 'solver 9', 'order.txt' u 1:4 ls 3 t 'solver 11', 'order.txt' u 1:5 ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Order Vs. Execution Time'
set ylabel 'Time [min]'
set output 'order_time.pdf'

plot 'order.txt' u 1:6 ls 1 t 'solver 8', 'order.txt' u 1:7 ls 2 t 'solver 9', 'order.txt' u 1:8 ls 3 t 'solver 11', 'order.txt' u 1:9 ls 4 t 'solver 12'

unset title
unset ylabel

set title 'Order Vs. Iterations'
set ylabel 'Iterations'
set output 'order_iter.pdf'

plot 'order.txt' u 1:10 ls 1 t 'solver 8', 'order.txt' u 1:11 ls 2 t 'solver 9', 'order.txt' u 1:12 ls 3 t 'solver 11', 'order.txt' u 1:13 ls 4 t 'solver 12'
