#!/bin/bash      

gnuplot -e "
set xrange [-20:450];
set yrange [-0:70];
set ticslevel 0;

set xlabel 'x';
set ylabel 'y';
set nokey;

filedata='out.txt';
n = system(sprintf('cat %s | wc -l', filedata));

do for [j=1:n-1] {
    set title 'time '.j;
    plot filedata u 2:3 every ::1::j w l lw 2, filedata u 2:3 every ::j::j w p pt 7 ps 2, 
         filedata u 4:5 every ::1::j w l lw 2, filedata u 4:5 every ::j::j w p pt 7 ps 2;
    pause 0.105;
};
"
