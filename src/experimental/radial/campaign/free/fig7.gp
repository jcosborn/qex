# Fig. 7 of arXiv:2510.03085: free overlap propagator G^(1,1)(t) vs the continuum.
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig7.gp
d = "output/radial/free/"
set terminal pngcairo size 900,650
set output d."fig7.png"
set logscale y
set xlabel "t"
set ylabel "G^{(1,1)}(t)"
set xrange [0:12]
set key top right
plot d."fig7_L1.tsv" u 1:2 w p pt 7 ps .5 t "L=1", \
     d."fig7_L2.tsv" u 1:2 w p pt 5 ps .5 t "L=2", \
     d."fig7_L4.tsv" u 1:2 w p pt 9 ps .5 t "L=4", \
     d."fig7_L8.tsv" u 1:2 w p pt 11 ps .5 t "L=8", \
     d."fig7_L1.tsv" u 1:3 w l lw 2 lc "black" t "(V.3) periodic images", \
     d."fig7_L1.tsv" u 1:4 w l dt 2 lc "gray" t "1/(16 pi sinh^2(t/2))"
