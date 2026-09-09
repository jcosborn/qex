# Fig. 12: gauge O(a^2) scaling, both projections of the (V.7) plane fit.
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig12.gp
d = "output/radial/free/"
set terminal pngcairo size 1200,520
set output d."fig12.png"
set multiplot layout 1,2
set key bottom right
set xlabel "abar_s^2"
set ylabel "Delta_0 - c_t a_t^2"
plot d."fig12_gauge_scaling.tsv" u 3:7 w p pt 7 t "lattice", \
     sqrt(2) w l dt 2 lc "gray" t "exact sqrt(2)"
set xlabel "a_t^2"
set ylabel "Delta_0 - c_s abar_s^2"
plot d."fig12_gauge_scaling.tsv" u 4:8 w p pt 7 t "lattice", \
     sqrt(2) w l dt 2 lc "gray" t "exact sqrt(2)"
unset multiplot
