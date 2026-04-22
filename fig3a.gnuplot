# SPDX-License-Identifier: MIT
# author: Emmanuel Le Trong <emmanuel.le-trong@cnrs-orleans.fr>
set logscale xy
set terminal pngcairo size 1920,1080
set key right
set format "%h"
set style data lines
set xlabel "R_b [m]"
set ylabel "F_{attach} / F_{detach} []"
set xrange [1e-7:1e-1]
set yrange [1e-6:1e9]
set arrow from 1e-7,1 to 1e-1,1 nohead dt 2

set output "fig3a_Asteroid.png"
set title "Asteroid"
plot  "fig_3a_Asteroid.csv" index 0 title "R_m = 4 μm" \
    , "" index 1 title "R_m = 16 μm" \
    , "" index 2 title "R_m = 64 μm" \
    , "" index 3 title "R_m = 256 μm" \
    , "" index 4 title "R_m = 1000 μm" \
    , "" index 5 title "R_m = 5000 μm" \

set output "fig3a_Vesta.png"
set title "Vesta"
plot  "fig_3a_Vesta.csv" index 0 title "R_m = 4 μm" \
    , "" index 1 title "R_m = 16 μm" \
    , "" index 2 title "R_m = 64 μm" \
    , "" index 3 title "R_m = 256 μm" \
    , "" index 4 title "R_m = 1000 μm" \
    , "" index 5 title "R_m = 5000 μm" \

set output "fig3a_Moon.png"
set title "Moon"
plot  "fig_3a_Moon.csv" index 0 title "R_m = 4 μm" \
    , "" index 1 title "R_m = 16 μm" \
    , "" index 2 title "R_m = 64 μm" \
    , "" index 3 title "R_m = 256 μm" \
    , "" index 4 title "R_m = 1000 μm" \
    , "" index 5 title "R_m = 5000 μm" \

set output "fig3a_Mars.png"
set title "Mars"
plot  "fig_3a_Mars.csv" index 0 title "R_m = 4 μm" \
    , "" index 1 title "R_m = 16 μm" \
    , "" index 2 title "R_m = 64 μm" \
    , "" index 3 title "R_m = 256 μm" \
    , "" index 4 title "R_m = 1000 μm" \
    , "" index 5 title "R_m = 5000 μm" \

set output "fig3a_Earth.png"
set title "Earth"
plot  "fig_3a_Earth.csv" index 0 title "R_m = 4 μm" \
    , "" index 1 title "R_m = 16 μm" \
    , "" index 2 title "R_m = 64 μm" \
    , "" index 3 title "R_m = 256 μm" \
    , "" index 4 title "R_m = 1000 μm" \
    , "" index 5 title "R_m = 5000 μm" \

