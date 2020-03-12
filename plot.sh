#!/bin/bash      

gnuplot -e "
set zrange [0:50];
set xrange [-150:400];
set yrange [-150:200];
set ticslevel 0;

set xlabel 'x';
set ylabel 'y';
set zlabel 'z';
col = 3;

filedata = 'out.txt';
n = system(sprintf('cat %s | wc -l', filedata));

do for [j=1:n-1] {
    set title 'time '.j;
    splot filedata u 2:3:4 every ::1::j w l lw 2, filedata u 2:3:4 every ::j::j w p pt 7 ps 2,
	  filedata u 5:6:7 every ::1::j w l lw 2, filedata u 5:6:7 every ::j::j w p pt 7 ps 2;
};
"
