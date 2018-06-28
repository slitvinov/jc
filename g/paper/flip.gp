sc = 0.75
set term pdf fontscale 0.80
set output "flip.pdf"

set key off

set xlabel "time"
set ylabel "eigenvalue"

set ytics 10

plot "flip" u 1:11 w dots lw 5 lc 0

set output