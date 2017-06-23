
set pm3d
set hidden3d

set view equal

splot 'ml_mos2_mlmc_L1_E_M_xx.bin' binary format="%double" array=(40,40) w l
 
set view 0,180
replot
