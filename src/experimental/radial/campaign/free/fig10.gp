# Fig. 10: T-symmetry: the overlap propagator folds about T/2, raw Wilson does not.
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig10.gp
d = "output/radial/free/"
set terminal pngcairo size 900,650
set output d."fig10.png"
set logscale y
set xlabel "t"
set ylabel "|G^{(1,1)}(t)|"
set xrange [0:12]
set key top center
plot d."fig10_L1.tsv" u 1:(abs($2)) w lp pt 7 ps .5 t "overlap", \
     d."fig10_L1.tsv" u 1:(abs($3)) w lp pt 5 ps .5 t "Wilson (raw D_W)"
