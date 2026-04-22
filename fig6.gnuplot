# SPDX-License-Identifier: MIT
# author: Emmanuel Le Trong <emmanuel.le-trong@cnrs-orleans.fr>
set logscale xy
set terminal pngcairo size 1920,1080
set key outside right
set format "%h"
set style data lines
set output "fig6.png"
set xlabel "R_m [m]"
set ylabel "R_b [m]"
set xrange [1e-7:1e0]
set yrange [1e-7:1e0]
plot  "fig_6_Asteroid_zero_buoyancy_b.csv" t "Asteroid" lc 2\
    , "fig_6_Asteroid_detachment1_b.csv" notitle lc 2\
    , "fig_6_Asteroid_detachment2_m.csv" notitle lc 2\
    , "fig_6_Vesta_zero_buoyancy_b.csv" t "Vesta" lc 1 \
    , "fig_6_Vesta_detachment1_b.csv" notitle lc 1 \
    , "fig_6_Vesta_detachment2_m.csv" notitle lc 1 \
    , "fig_6_Moon_zero_buoyancy_b.csv" t "Moon" lc 3 \
    , "fig_6_Moon_detachment1_b.csv" notitle lc 3 \
    , "fig_6_Moon_detachment2_m.csv" notitle lc 3 \
    , "fig_6_Mars_zero_buoyancy_b.csv" t "Mars" lc 4 \
    , "fig_6_Mars_detachment1_b.csv" notitle lc 4 \
    , "fig_6_Mars_detachment2_m.csv" notitle lc 4 \
    , "fig_6_Earth_zero_buoyancy_b.csv" t "Earth" lc 5 \
    , "fig_6_Earth_detachment1_b.csv" notitle lc 5 \
    , "fig_6_Earth_detachment2_m.csv" notitle lc 5 \
