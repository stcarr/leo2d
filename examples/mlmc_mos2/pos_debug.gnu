set style circle radius .001

x_c = -1;
y_c = -1;
axis_lim = 15;

set xrange [-axis_lim+x_c:axis_lim+x_c]
set yrange [-axis_lim+y_c:axis_lim+y_c]
set zrange [-3:3]

set view 90,0

tar_orb = 231

splot 'ml_mos2_mlmc_pos.dat' u 2:3:4 with points, \
	'ml_mos2_mlmc_inter_pos.dat' using 1:2:($7 ==tar_orb ? $3 : 1/0):($4-$1):($5-$2):($6-$3) with vectors nohead, \
	'ml_mos2_mlmc_intra_pos.dat' using 1:2:($7 ==tar_orb ? $3 : 1/0):($4-$1):($5-$2):($6-$3) with vectors nohead

