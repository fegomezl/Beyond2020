set g
set term pdf
set key opaque
set key t r

set ls 1 lc rgb 'blue' lt 1 pt 1 lw 1 ps 1
set ls 2 lc rgb 'red' lt 2 pt 2 lw 1 ps 1
set ls 3 lc rgb 'black' lt 6 pt 6 lw 1 ps 1

set title 'Total mean error vs Border width'
set xl 'Border (u.a)'
set yl 'Error (u.a)'

set xr [:]
set yr [:]
set log

set o "border_step3.pdf"

plot 'border.txt' u 2:3 t 'W error' ls 1, \
     'border.txt' u 2:4 t '{/Symbol y} error' ls 2, 

set o "inv_r.pdf"

set key t l
set title 'Total mean error vs inverse r constant'
set xl 'Inv_r (u.a)'

set xr [:]
set yr [:]


plot 'inv_r.txt' u 2:3 t 'W error' ls 1, \
     'inv_r.txt' u 2:4 t '{/Symbol y} error' ls 2

set o "refs_step3.pdf"

set key t l
set title 'Total mean error vs Refinements'
set xl 'Mesh Size (u.a)'

plot 'refs.txt' u 2:3 t 'W error' ls 1, \
     'refs.txt' u 2:4 t '{/Symbol y} error' ls 2

set o "order_step3.pdf"

set key t l
set title 'Total mean error vs Order'
unset logscale x
unset xl
set xl 'Order'
set xtics 1

plot 'order.txt' u 1:3 t 'W error' ls 1, \
     'order.txt' u 1:4 t '{/Symbol y} error' ls 2
