set style circle radius .001

filename = 'mlgrphn'

x_c = 2;
y_c = 5;
axis_lim = 6;

set xrange [-axis_lim+x_c:axis_lim+x_c]
set yrange [-axis_lim+y_c:axis_lim+y_c]
set zrange [-1:5]

set view 90,0

tar_orb = 0

splot filename.'_pos.dat' u 2:3:4 with points, \
	filename.'_inter_pos.dat' using 1:2:($7 ==tar_orb ? $3 : 1/0):($4-$1):($5-$2):($6-$3) with vectors nohead, \
	filename.'_intra_pos.dat' using 1:2:($7 ==tar_orb ? $3 : 1/0):($4-$1):($5-$2):($6-$3) with vectors nohead, \
	'supercell.dat' using 1:2:3:($4-$1):($5-$2):($6-$3) with vectors nohead

set view 0,0
replot

#pause -1
