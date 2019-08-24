
# first run transpose.sh on *.cheb file!!
# sudo /bin/bash /home/stephen/devspace/harvard/LEO2D_devbranch/realspace/transpose.sh > mlgrphn_transpose.cheb

num_samps = 1;

plot 'mlgrphn_transpose.cheb' u (sum [col=1:2*num_samps] column(col)) w lines

