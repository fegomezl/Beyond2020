set g
set title ''
set xlabel 'Time (min)'
set ylabel 'Interface (mm)'

set ls 1 lc rgb "black" pt 7 ps 0.2
set ls 2 lc rgb "red" lt 1 lw 1

file = 'results/graph/data.txt'

set datafile separator ","
set key autotitle columnhead
set term pdf
set o 'results/regression.pdf'
set fit quiet

f(x) = A*x**B
A=1
B=1

fit f(x) file every ::1 u 1:2:3 yerrors via A,B

Fit = sprintf(" {/:Bold Parameters} \n y = ax^b \n a = %g +/- %g \n b = %g +/- %g", A, A_err, B, B_err)

set obj 2 rect from graph 0, 1 to graph 0.50, 0.77 fc rgb "white"
set lab 2 Fit at graph 0, 0.97


plot file every ::1 u 1:2:3 w yerrorbars ls 1 t 'data', f(x) w l t 'regression' ls 2
