# Fig. 11: free gauge current correlator G_g(t) vs the continuum (V.14).
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig11.gp
d = "output/radial/free/"
set terminal pngcairo size 900,650
set output d."fig11.png"
set logscale y
set xlabel "t"
set ylabel "G_g(t)"
set xrange [0:12]
set key top right
plot d."fig11_L1.tsv" u 1:2 w p pt 7 ps .5 t "L=1", \
     d."fig11_L2.tsv" u 1:2 w p pt 5 ps .5 t "L=2", \
     d."fig11_L4.tsv" u 1:2 w p pt 9 ps .5 t "L=4", \
     d."fig11_L8.tsv" u 1:2 w p pt 11 ps .5 t "L=8", \
     d."fig11_L1.tsv" u 1:3 w l lw 2 lc "black" t "(V.14) periodic images"
