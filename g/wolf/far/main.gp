set term pdf
set output "main.pdf"

unset key
unset tics
unset border
set size ratio -1
rot_x = 84; rot_y = 40; sc = 1.35
set  view rot_x, rot_y, sc, sc

splot "q" w l lw 2 lc  rgb "#00dddddd", \
      "v" w l lw 3 lc  "black", \
      "p"  w p pt 7 ps 2 lc 0
set output