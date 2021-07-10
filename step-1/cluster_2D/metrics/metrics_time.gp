set g

set title 'CPU time vs h_{min} for order 1'
set xr [:]
set yr [:]
set xl 'h_{min}' font "arial,14"
set yl 'time [s]' font "arial, 14"

set term pdf font "arial, 14"

f1(x) = A1*x**B1
f2(x) = A2*x**B2
f3(x) = A3*x**B3

fit f1(x) 'metrics/2D_structured_1.txt' u 1:2 via A1, B1
fit f2(x) 'metrics/2D_unstructured_1.txt' u 1:2 via A2, B2
fit f3(x) '../cluster_3D/metrics/3D_1.txt' u 1:2 via A3, B3

set o "metrics/time_step1-order1.pdf"
plot 'metrics/2D_structured_1.txt' u 1:2 t '2D structured' pt 1 lc rgb "blue",\
'metrics/2D_unstructured_1.txt' u 1:2 t '2D unstructured' pt 2 lc rgb "red",\
'../cluster_3D/metrics/3D_1.txt' u 1:2 t '3D structured' pt 6 lc rgb "black"

set title 'time vs h_{min} for order 1 (log scale)'
set log
set xr []
set yr [:]
set o "metrics/time_step1-order1_log.pdf"

rep
