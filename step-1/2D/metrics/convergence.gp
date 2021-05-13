set g

set title 'error vs h_{min} for order 1'
set xr [0:0.85]
set yr [0:5e-5]
set xl 'h_{min}'
set yl 'error'

set term pdf

f1(x) = A1*x**B1
f2(x) = A2*x**B2
f3(x) = A3*x**B3

fit f1(x) 'metrics/convergence-structured.txt' u 2:3 via A1, B1
fit f2(x) 'metrics/convergence-unstructured.txt' u 2:3 via A2, B2
fit f3(x) 'metrics/convergence-3D.txt' u 2:3 via A3, B3

set o "metrics/convergence_step1-order1.pdf"
plot 'metrics/convergence-structured.txt' u 2:3 t '2D structured' pt 1 lc rgb "blue", f1(x) lc rgb "blue" notitle,\
'metrics/convergence-unstructured.txt' u 2:3 t '2D unstructured' pt 2 lc rgb "red", f2(x) lc "red" notitle,\
'metrics/convergence-3D.txt' u 2:3 t '3D' pt 6 lc rgb "black", f3(x) lc rgb "black" notitle

set log
set xr [1e-3:1]
set yr [1e-11:1e-3]
set o "metrics/convergence_step1-order1_log.pdf"

rep
