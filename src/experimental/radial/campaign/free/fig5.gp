# Fig. 5: the flat-limit spectrum (IV.8) over the curved one (same panels as Fig. 4).
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig5.gp
d = "output/radial/free/"
set terminal pngcairo size 1400,900
set output d."fig5.png"
set multiplot layout 2,3
set xlabel "Re"
set ylabel "Im"
set title "L=1, Lt=24"
plot d."fig5_flat_L1_Lt24.tsv" u 1:2 w p pt 7 ps .2 lc "gray" t "flat (IV.8)", \
     d."fig4_L1_Lt24.tsv" u 1:2 w p pt 7 ps .3 lc "red" t "curved"
set title "L=2, Lt=24"
plot d."fig5_flat_L2_Lt24.tsv" u 1:2 w p pt 7 ps .2 lc "gray" not, \
     d."fig4_L2_Lt24.tsv" u 1:2 w p pt 7 ps .3 lc "red" not
set title "L=4, Lt=24"
plot d."fig5_flat_L4_Lt24.tsv" u 1:2 w p pt 7 ps .2 lc "gray" not, \
     d."fig4_L4_Lt24.tsv" u 1:2 w p pt 7 ps .3 lc "red" not
set title "L=2, Lt=16"
plot d."fig5_flat_L2_Lt16.tsv" u 1:2 w p pt 7 ps .2 lc "gray" not, \
     d."fig4_L2_Lt16.tsv" u 1:2 w p pt 7 ps .3 lc "red" not
set title "L=2, Lt=48"
plot d."fig5_flat_L2_Lt48.tsv" u 1:2 w p pt 7 ps .2 lc "gray" not, \
     d."fig4_L2_Lt48.tsv" u 1:2 w p pt 7 ps .3 lc "red" not
unset multiplot
