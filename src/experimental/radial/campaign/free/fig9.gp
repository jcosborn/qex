# Fig. 9: fermion O(a^2) scaling, both projections of the (V.7) plane fit.
# proj_s = Delta_0 - c_t a_t^2 (vs abar^2), proj_t = Delta_0 - c_s abar^2 (vs a_t^2).
# Run from the repo root:  gnuplot src/experimental/radial/campaign/free/fig9.gp
d = "output/radial/free/"
set terminal pngcairo size 1200,520
set output d."fig9.png"
set multiplot layout 1,2
set key top left
set xlabel "abar_s^2"
set ylabel "Delta_0 - c_t a_t^2"
plot d."fig9_fermion_scaling.tsv" u 3:7 w p pt 7 t "lattice", \
     1 w l dt 2 lc "gray" t "exact 1"
set xlabel "a_t^2"
set ylabel "Delta_0 - c_s abar_s^2"
plot d."fig9_fermion_scaling.tsv" u 4:8 w p pt 7 t "lattice", \
     1 w l dt 2 lc "gray" t "exact 1"
unset multiplot
