#set terminal png
#set output "plot.png"

set terminal pdf size 6,2.25
set output "plot_noTT_long.pdf"
#set terminal x11
set multiplot layout 1,2

unset offsets
unset key
set offsets graph 0.07, graph 0.07, 0, 0

set xlabel "Timesteps"
set ylabel "Mean Concentration"
set ylabel offset 1,0

set xrange [0:200]
set yrange [-0.2:4]

set format y "%.1f"
set grid ytics lc rgb "#bbbbbb" lw 1 lt 0
set grid xtics lc rgb "#bbbbbb" lw 1 lt 0

LMARGIN = "set lmargin at screen 0.06; set rmargin at screen 0.515"
RMARGIN = "set lmargin at screen 0.515; set rmargin at screen 0.99"

@LMARGIN
set title "Tensor-Train Format"
plot for [i=1:5] "mean_TT.dat" u 0:i w l t "protein ".i
#set output

@RMARGIN
set title "Simple Tensor Format"
set ytics scale 0
set format y ""
set grid
set offsets graph 0.07, graph 0.07, 0, 0
unset ylabel
set key right top
plot for [i=1:5] "mean_noTT.dat" u 0:i w l t "protein ".i

unset multiplot
unset output
pause -1 "Hit any key to continue"