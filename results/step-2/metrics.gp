set g
set term pdf
set key opaque
set key t r

set ls 1 lc rgb 'blue' lt 1 pt 1 lw 1 ps 1
set ls 2 lc rgb 'red' lt 2 pt 2 lw 1 ps 1
set ls 3 lc rgb 'black' lt 6 pt 6 lw 1 ps 1

set title 'CPU time vs Mesh size for order 1 (log-log)'
set xl 'Mesh size (u.a)'
set yl 'Time (s)'

set xr [:]
set yr [:]
set log

set o "time_step1.pdf"

plot '2D_structured_1.txt' u 1:2 t '2D structured' ls 1, \
     '2D_unstructured_1.txt' u 1:2 t '2D unstructured' ls 2, \
     '3D_1.txt' u 1:2 t '3D' ls 3

set o "error_step1.pdf"

set key t l
set title 'Error vs Mesh size for order 1 (log-log)'
set xl 'Mesh size (u.a)'
set yl 'Error (u.a)'

set xr [:]
set yr [:]
set log

plot '2D_structured_1.txt' u 1:3 t '2D structured' ls 1, \
     '2D_unstructured_1.txt' u 1:3 t '2D unstructured' ls 2, \
     '3D_1.txt' u 1:3 t '3D' ls 3
