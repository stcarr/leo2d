set style circle radius .001

filename = 'mlgrphn'

x_c = -1;
y_c = -1;
axis_lim = 50;

set xrange [-axis_lim+x_c:axis_lim+x_c]
set yrange [-axis_lim+y_c:axis_lim+y_c]
set zrange [-1:5]

set view 90,0

tar_orb = 0

splot filename.'_pos.dat' u 2:3:4 with points, \
	filename.'_inter_pos.dat' using 1:2:($7 ==tar_orb ? $3 : 1/0):($4-$1):($5-$2):($6-$3) with vectors nohead, \
	filename.'_intra_pos.dat' using 1:2:($7 ==tar_orb ? $3 : 1/0):($4-$1):($5-$2):($6-$3) with vectors nohead

pause -1
