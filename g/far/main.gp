set term pdf
set output "main.pdf"

unset key
unset tics
unset border
set size ratio -1
rot_x = 70; rot_y = 110; sc = 1.33
set  view rot_x, rot_y, sc, sc

splot "q"  w l lw 2 lc  rgb "#00dddddd", \
      "q0" w l lw 4 lc "red"
set output