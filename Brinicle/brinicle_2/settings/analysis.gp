set g
set title ''
set xlabel 'Iteration'
set ylabel 'Dt'

set ls 1 lc rgb "black" pt 7 ps 0.2
set ls 2 lc rgb "red" lt 1 lw 1

file = 'results/progress.txt'

set key autotitle columnhead
set key r b
set term pdf
set o 'results/analysis.pdf'

plot file u 1:2 ls 1 t 'Dt'
