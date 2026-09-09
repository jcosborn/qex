# Fig. 4: free D_W spectra (raw D_lat units), T = 4: L scan at Lt = 24 and Lt scan
# at L = 2.  Run from the repo root:
#   gnuplot src/experimental/radial/campaign/free/fig4.gp
d = "output/radial/free/"
set terminal pngcairo size 1400,900
set output d."fig4.png"
set multiplot layout 2,3
set xlabel "Re"
set ylabel "Im"
set title "L=1, Lt=24"
plot d."fig4_L1_Lt24.tsv" u 1:2 w p pt 7 ps .3 not
set title "L=2, Lt=24"
plot d."fig4_L2_Lt24.tsv" u 1:2 w p pt 7 ps .3 not
set title "L=4, Lt=24"
plot d."fig4_L4_Lt24.tsv" u 1:2 w p pt 7 ps .3 not
set title "L=2, Lt=16"
plot d."fig4_L2_Lt16.tsv" u 1:2 w p pt 7 ps .3 not
set title "L=2, Lt=48"
plot d."fig4_L2_Lt48.tsv" u 1:2 w p pt 7 ps .3 not
unset multiplot
