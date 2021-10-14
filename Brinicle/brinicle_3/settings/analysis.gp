set g
set title 'Dt vs iterations'
set xlabel 'Iteration'
set ylabel 'Dt(10^{(-n)})'

set ls 1 lc rgb "black" lt 1 lw 1 pt 2 ps 0.5

file = 'results/progress.txt'

set key autotitle columnhead
set key r b
set term pdf
set o 'results/analysis.pdf'

plot[0:][0:] file u 1:(-log10($2)) w l ls 1 t 'Dt'
