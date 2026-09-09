# Fig. 8: fermion Delta_eff(t), T = 16, Lt = 168.
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig8.gp
d = "output/radial/free/"
set terminal pngcairo size 900,650
set output d."fig8.png"
set xlabel "t"
set ylabel "Delta_{eff}(t)"
set xrange [0:8]
set yrange [0.8:1.6]
set key bottom right
plot d."fig8_L1.tsv" u 1:2 w lp pt 7 ps .4 t "L=1", \
     d."fig8_L2.tsv" u 1:2 w lp pt 5 ps .4 t "L=2", \
     d."fig8_L4.tsv" u 1:2 w lp pt 9 ps .4 t "L=4", \
     d."fig8_L8.tsv" u 1:2 w lp pt 11 ps .4 t "L=8", \
     d."fig8_L1.tsv" u 1:3 w l lw 2 lc "black" t "continuum (V.3)", \
     1 w l dt 2 lc "gray" t "Delta_0 = 1"
