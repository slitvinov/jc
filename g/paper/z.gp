sc = 0.75
set term pdf fontscale 0.80
set output "z.pdf"

set key off

set xlabel "time"
set ylabel "z"

set ytics 0.25

plot "q" u 1:4 w l lw 5 lc 0

set output