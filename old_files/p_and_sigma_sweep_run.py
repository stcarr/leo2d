#Script for running multiple instances of LEO2D on Harvard Odyssey

import os
import sys
import string
import numpy as np
import math

#main program starts here

main_dir=os.getcwd()+'/'

sigma_list = [0.32, 0.16, 0.08, 0.04, 0.02, 0.01]

angle = 3

dir_name_mom = main_dir+'mom_v2'
os.mkdir(dir_name_mom)
os.chdir(dir_name_mom)

for sigma in sigma_list:
	
	R = 7 # int(math.ceil(1.0*math.log(1/sigma)))
	N = int(math.ceil((1.0/3.0)*(1/sigma)*math.log(1/sigma)))
	print(sigma)
	print(R)
	print(N)
	dir_name=main_dir+string.replace("%.04lf" % angle,'.','p')+'/'
	os.system('cp ../hstruct.in hstruct_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'.in')
	#os.system('cp ../run_script.sh .')
	os.system('sed -i \'s/%NAME/'+'blg_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'/g\' hstruct_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'.in')
	os.system('sed -i \'s/%SIZE/'+str(R)+'/g\' hstruct_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'.in')
	os.system('sed -i \'s/%N/'+str(N)+'/g\' hstruct_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'.in')
	os.system('sed -i \'s/%ANG/'+str(angle)+'/g\' hstruct_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'.in')
	#os.system('sbatch run_script.sh')   
	os.system('mpirun -n 2 LEO2D hstruct_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'.in | tee hstruct_sigma_'+string.replace("%0.4lf" % sigma,'.','p')+'.out')

