
set pm3d
set hidden3d

set view equal

splot 'L1_MXX_NO_K.bin' binary format="%double" array=(40,40) w l
 
set view 0,180
replot
