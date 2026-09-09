# Fig. 6: D_W vs D_ov (T = 4, L = 4, Lt = 24, M = 1), raw spectra and the
# generalized eigenvalues with the corrected volume weight diag(volbar/volw).
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig6.gp
d = "output/radial/free/"
set terminal pngcairo size 1200,900
set output d."fig6.png"
set multiplot layout 2,2
set xlabel "Re"
set ylabel "Im"
set title "raw D_W"
plot d."fig6_dw.tsv" u 1:2 w p pt 7 ps .3 not
set title "raw D_ov (Ginsparg-Wilson circle)"
set object 1 circle at 1,0 size 1 fs empty border lc "gray"
plot d."fig6_dov.tsv" u 1:2 w p pt 7 ps .3 not
unset object 1
set title "gen. eig. of D_W (weight volbar/volw)"
plot d."fig6_dw_gen.tsv" u 1:2 w p pt 7 ps .3 not
set title "gen. eig. of D_ov (weight volbar/volw)"
plot d."fig6_dov_gen.tsv" u 1:2 w p pt 7 ps .3 not
unset multiplot
