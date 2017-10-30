import csv, os, sys
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

######################
# bold_confounds = [.tsv file containing all bold confounds $OAK/projects/[sunetID]/fmriprep/sub-[PIDvs]/func/*.tsv 
######################

def main():
	if len(sys.argv) == 3:

		script = sys.argv[0]
                tsvfile_sub = sys.argv[1]
		tsvfile_task = sys.argv[2]
		tsvfile = os.getcwd()
		tsvfile = tsvfile.split('/')
		tsvfile = '/'.join(tsvfile[:-2])+'/fmriprep/sub-%s/func/sub-%s_task-%s_bold_confounds.tsv'%(tsvfile_sub,tsvfile_sub,tsvfile_task)
                get_bold_conf(tsvfile)
	
	else:
		print("Use:"+"\n"+"\n"+"python plt_bold_confounds.py [subject PID, visit, session (e.g. '354411')] [task (e.g. 'resting_state_1']"+"\n")

def get_bold_conf(tsvfile):
	
	pwd = os.getcwd()
	bold_confound_path = tsvfile
	tsvfile = tsvfile.split('-')

	df = pd.read_csv(bold_confound_path,sep='\t')

	# Last six columns of the tsv are the movement confounds 

	x = df['X']
	y = df['Y']
	z = df['Z']
	rotx = df['RotX']
	roty = df['RotY']
	rotz = df['RotZ']

	# Plot x y z

	xplt = plt.plot(x,'r')
	yplt = plt.plot(y,'y')
	zplt = plt.plot(z,'b')
	plt.ylabel('Motion (mm)')
	plt.xlabel('Timesteps (sub-bricks)')
	plt.legend((xplt[0],yplt[0],zplt[0]), ('X motion','Y motion','Z motion'))
	plt.savefig(pwd+'/'+'xyzmotion.png')
	plt.show()

	# Plot rotations (x y z)

	xrot = plt.plot(rotx,'c')
	yrot = plt.plot(roty,'g')
	zrot = plt.plot(rotz,'m')
	plt.ylabel('Motion (mm)')
	plt.xlabel('Timesteps (sub-bricks)')
	plt.legend((xrot[0],yrot[0],zrot[0]), ('X rotation motion','Y rotation motion','Z rotation motion'))
	plt.savefig(pwd+'/'+'rotational_motion.png')
	plt.show()

if __name__ == '__main__':
	main()
