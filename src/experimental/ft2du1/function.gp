# Usage:
#   gnuplot -e "files='sine.dat bspline.dat'; names='sine bspline'; out='maps.png'" function.gp
# Select one data set from each file with, for example, indices='0 2'.

if (!exists("files")) files = "map-function.dat"
if (!exists("names")) names = files
if (!exists("indices")) indices = ""
if (!exists("out")) out = "map-functions.png"
if (!exists("titleText")) titleText = "Analytic scalar map functions"

nfile = words(files)
idx(i) = words(indices) >= i ? int(word(indices, i)) : 0

set terminal pngcairo size 1800,1350 enhanced font "Helvetica,15"
set output out
set datafile commentschars "#"
set encoding utf8
set border lw 1.2
set xrange [-pi:pi]
set xtics ("-pi" -pi, "-pi/2" -pi/2, "0" 0, "pi/2" pi/2, "pi" pi)
set xlabel "auxiliary scalar input x"
set tics out
set grid xtics ytics lc rgb "#d8d8d8" lw 0.7
set key inside top left vertical opaque box lc rgb "#b8b8b8" samplen 2

set linetype 1 lc rgb "#0072B2" lw 2.4
set linetype 2 lc rgb "#D55E00" lw 2.4
set linetype 3 lc rgb "#009E73" lw 2.4
set linetype 4 lc rgb "#CC79A7" lw 2.4
set linetype 5 lc rgb "#E69F00" lw 2.4
set linetype 6 lc rgb "#56B4E9" lw 2.4
set linetype 7 lc rgb "#8C564B" lw 2.4
set linetype 8 lc rgb "#6F4E7C" lw 2.4

set multiplot layout 2,2 rowsfirst title titleText font ",19" margins 0.07,0.97,0.075,0.91 spacing 0.09,0.12

set title "Forward function"
set ylabel "g(x)"
plot x with lines lc rgb "#777777" dt 2 lw 1.4 title "identity", \
  for [i=1:nfile] word(files,i) index idx(i) using 4:5 with lines ls i title word(names,i)

set title "Periodic displacement"
set ylabel "g(x)-x"
plot 0 with lines lc rgb "#777777" dt 2 lw 1.4 notitle, \
  for [i=1:nfile] word(files,i) index idx(i) using 4:6 with lines ls i title word(names,i)

set title "Scalar Jacobian"
set ylabel "g'(x)"
plot 1 with lines lc rgb "#777777" dt 2 lw 1.4 notitle, \
  for [i=1:nfile] word(files,i) index idx(i) using 4:7 with lines ls i title word(names,i)

set title "Jacobian force"
set ylabel "d log(g'(x)) / dx"
plot 0 with lines lc rgb "#777777" dt 2 lw 1.4 notitle, \
  for [i=1:nfile] word(files,i) index idx(i) using 4:10 with lines ls i title word(names,i)

unset multiplot
unset output
