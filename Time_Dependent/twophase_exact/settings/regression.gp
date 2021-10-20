set g
set title 'Interface position vs time'
set xlabel 'Time (min)'
set ylabel 'Interface (mm)'

set ls 1 lc rgb "black" pt 7 ps 0.2
set ls 2 lc rgb "red" lt 1 lw 1

file = 'results/data.txt'

set datafile separator ","
set key autotitle columnhead
set key r b
set term pdf
set o 'results/regression.pdf'
set fit quiet

f(x) = A*x**B
A=2.55503
B=0.5

fit[0:] f(x) file u 1:2 via A,B

Fit = sprintf(" {/:Bold Parameters} \n y = ax^b \n a = %.5f ± %.5f \n b = %.5f ± %.5f \n Δy_{max} = 0.005 mm", A, A_err, B, B_err)

set obj 2 rect from graph 0, 1 to graph 0.3, 0.75 fc rgb "white"
set lab 2 Fit at graph 0, 0.97


plot[0:1] file u 1:2 ls 1 t 'Interface', f(x) w l t 'Fit' ls 2
