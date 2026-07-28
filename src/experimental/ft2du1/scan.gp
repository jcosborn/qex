# Usage:
#   gnuplot -e "file='scan.dat'; out='scan.png'; geom='link2'; beta=6" scan.gp
# Optional common comparison ranges: potMax, logjMax, forceMax, stiffMax.
# Non-potential contours use each data range; set *ContMin/*ContMax to override.
# Set labelOther=0 to omit labels on the non-potential contours.

if (!exists("file")) file = "scan.dat"
if (!exists("out")) out = file . ".png"
if (!exists("geom")) geom = "plaq4"
if (!exists("beta")) beta = 6.0
if (!exists("name")) name = file
if (!exists("ncont")) ncont = 12
if (!exists("ncontOther")) ncontOther = 6
if (!exists("labelStart")) labelStart = 8
if (!exists("labelEvery")) labelEvery = 72
if (!exists("labelStartOther")) labelStartOther = 24
if (!exists("labelEveryOther")) labelEveryOther = 260
if (!exists("labelOther")) labelOther = 1

if (!exists("potMin")) potMin = 0.0
if (!exists("potMax")) potMax = 2.0*beta
if (!exists("logjMax")) logjMax = 0.5
if (!exists("forceMax")) forceMax = 4.0*beta
if (!exists("stiffMax")) stiffMax = 24.0*beta

stats file using 6 nooutput
logjDataMin = STATS_min
logjDataMax = STATS_max
stats file using (sqrt($10)) nooutput
forceDataMin = STATS_min
forceDataMax = STATS_max
stats file using 11 nooutput
stiffDataMin = STATS_min
stiffDataMax = STATS_max

if (!exists("logjContMin")) logjContMin = logjDataMin
if (!exists("logjContMax")) logjContMax = logjDataMax
if (!exists("forceContMin")) forceContMin = forceDataMin
if (!exists("forceContMax")) forceContMax = forceDataMax
if (!exists("stiffContMin")) stiffContMin = stiffDataMin
if (!exists("stiffContMax")) stiffContMax = stiffDataMax
logjLabelSpan = logjContMax-logjContMin
if (logjLabelSpan <= 0.0) logjLabelSpan = 1.0
logjLabelY(z) = -0.65*pi+1.3*pi*(z-logjContMin)/logjLabelSpan

if (geom eq "link2") {
  xlab = "P_{+}"
  ylab = "P_{-}"
} else {
  xlab = "P_c"
  ylab = "P_n  (P_0=P_1=P_2=P_3)"
}

set terminal pngcairo size 1800,1450 enhanced font "Helvetica,15"
set output out
set datafile commentschars "#"
set encoding utf8
set border lw 1.2
set key off
set size ratio -1
set xrange [-pi:pi]
set yrange [-pi:pi]
set xtics ("-pi" -pi, "-pi/2" -pi/2, "0" 0, "pi/2" pi/2, "pi" pi)
set ytics ("-pi" -pi, "-pi/2" -pi/2, "0" 0, "pi/2" pi/2, "pi" pi)
set xlabel xlab
set ylabel ylab
set tics out
set style textbox 1 opaque margins 0.28,0.14 fc rgb "#303030" noborder

# Extract common-level contours before entering the multiplot.
set contour base
unset surface
set cntrparam levels incremental potMin,(potMax-potMin)/ncont,potMax
set table $potContours
splot file using 1:2:9
unset table

logjContStep = (logjContMax-logjContMin)/(ncontOther+1.0)
if (logjContStep <= 0.0) logjContStep = 2.0*logjMax/(ncontOther+1.0)
set cntrparam levels incremental logjContMin+logjContStep,logjContStep,logjContMax-logjContStep
set table $logjContours
splot file using 1:2:6
unset table

forceContStep = (forceContMax-forceContMin)/(ncontOther+1.0)
if (forceContStep <= 0.0) forceContStep = forceMax/(ncontOther+1.0)
set cntrparam levels incremental forceContMin+forceContStep,forceContStep,forceContMax-forceContStep
set table $forceContours
splot file using 1:2:(sqrt($10))
unset table

stiffContStep = (stiffContMax-stiffContMin)/(ncontOther+1.0)
if (stiffContStep <= 0.0) stiffContStep = stiffMax/(ncontOther+1.0)
set cntrparam levels incremental stiffContMin+stiffContStep,stiffContStep,stiffContMax-stiffContStep
set table $stiffContours
splot file using 1:2:11
unset table

unset contour
set surface

set multiplot layout 2,2 rowsfirst title name font ",19" margins 0.055,0.965,0.07,0.925 spacing 0.085,0.10

set title "Effective potential per plaquette"
set cblabel "Delta S_eff / N_plaq"
set cbrange [potMin:potMax]
set palette defined (0 "#440154", 0.25 "#3b528b", 0.50 "#21918c", 0.75 "#5ec962", 1 "#fde725")
plot file using 1:2:9 with image, \
  $potContours using 1:2:3 with lines lc rgb "#f4f4f4" lw 0.8, \
  $potContours every labelEvery::labelStart using 1:2:(sprintf("%.1f",$3)) \
    with labels font "Helvetica,9" tc rgb "#ffffff" boxed bs 1

set title "Local log-Jacobian per active variable"
set cblabel "log(det J) / N_var"
set cbrange [-logjMax:logjMax]
set palette defined (0 "#3b4cc0", 0.25 "#9ebeff", 0.50 "#f7f7f7", 0.75 "#f6a385", 1 "#b40426")
if (labelOther) {
  if (geom eq "plaq4") {
    plot file using 1:2:6 with image, \
      $logjContours using 1:2:3 with lines lc rgb "#303030" lw 0.7, \
      $logjContours every labelEveryOther::labelStartOther using 1:(logjLabelY($3)):(sprintf("%.2g",$3)) \
        with labels font "Helvetica,9" tc rgb "#ffffff" boxed bs 1
  } else {
    plot file using 1:2:6 with image, \
      $logjContours using 1:2:3 with lines lc rgb "#303030" lw 0.7, \
      $logjContours every labelEveryOther::labelStartOther using 1:2:(sprintf("%.2g",$3)) \
        with labels font "Helvetica,9" tc rgb "#ffffff" boxed bs 1
  }
} else {
  plot file using 1:2:6 with image, \
    $logjContours using 1:2:3 with lines lc rgb "#303030" lw 0.7
}

set title "HMC force magnitude"
set cblabel "sqrt(force2)"
set cbrange [0:forceMax]
set palette defined (0 "#000004", 0.20 "#3b0f70", 0.45 "#8c2981", 0.70 "#de4968", 0.87 "#fe9f6d", 1 "#fcfdbf")
if (labelOther) {
  plot file using 1:2:(sqrt($10)) with image, \
    $forceContours using 1:2:3 with lines lc rgb "#f4f4f4" lw 0.7, \
    $forceContours every labelEveryOther::labelStartOther using 1:2:(sprintf("%.1f",$3)) \
      with labels font "Helvetica,9" tc rgb "#ffffff" boxed bs 1
} else {
  plot file using 1:2:(sqrt($10)) with image, \
    $forceContours using 1:2:3 with lines lc rgb "#f4f4f4" lw 0.7
}

set title "Largest absolute force-gradient eigenvalue"
set cblabel "stiffMax"
set cbrange [0:stiffMax]
set palette defined (0 "#00204c", 0.20 "#31446b", 0.45 "#666970", 0.70 "#a38f60", 0.87 "#ddb64a", 1 "#ffea46")
if (labelOther) {
  plot file using 1:2:11 with image, \
    $stiffContours using 1:2:3 with lines lc rgb "#f4f4f4" lw 0.7, \
    $stiffContours every labelEveryOther::labelStartOther using 1:2:(sprintf("%.0f",$3)) \
      with labels font "Helvetica,9" tc rgb "#ffffff" boxed bs 1
} else {
  plot file using 1:2:11 with image, \
    $stiffContours using 1:2:3 with lines lc rgb "#f4f4f4" lw 0.7
}

unset multiplot
unset output
