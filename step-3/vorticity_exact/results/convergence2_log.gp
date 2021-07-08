set terminal jpeg enhanced
set termoption enhanced
set grid
set style line 1 lc rgb "blue" lt 1 lw 2 pt 7 ps 0.5
set style line 2 lc rgb "red" lt 1 lw 2 pt 7 ps 0.5
set style line 3 lc rgb "green" lt 1 lw 2 pt 7 ps 0.5
set style line 4 lc rgb "orange" lt 1 lw 2 pt 7 ps 0.5

set title 'Refinaments Vs. Total Mean Error Case 1'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set output 'refs_error_log0.jpg'
set logscale x
set logscale y

plot 'refs0.txt' u 2:3 w l ls 1 t 'Error for w', 'refs0.txt' u 2:4 w l ls 2 t 'Error for {/Symbol y}'

unset title
unset ylabel

set title 'Order Vs. Total Mean Error Case 1'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set output 'order_error_log0.jpg'

plot 'order0.txt' u 2:3 w l ls 1 t 'Error for w', 'order0.txt' u 2:4 w l ls 2 t 'Error for {/Symbol y}'

set title 'Refinaments Vs. Total Mean Error Case 2'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set output 'refs_error_log1.jpg'
set logscale x
set logscale y

plot 'refs1.txt' u 2:3 w l ls 1 t 'Error for w', 'refs1.txt' u 2:4 w l ls 2 t 'Error for {/Symbol y}'

unset title
unset ylabel


set title 'Order Vs. Total Mean Error Case 2'
set xlabel 'Size'
set ylabel 'Total Mean Error'
set output 'order_error_log1.jpg'

plot 'order1.txt' u 2:3 w l ls 1 t 'Error for w', 'order1.txt' u 2:4 w l ls 2 t 'Error for {/Symbol y}'

