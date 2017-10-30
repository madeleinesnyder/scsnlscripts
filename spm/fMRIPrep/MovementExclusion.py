# MovementExclusion.py

# Generate a .csv of analogous movement stats for a given batch of preprocessed functionals

import os, sys, csv
import pandas as pd
import numpy as np
import subprocess
from subprocess import call

def extract_stats(tsvfile):
        df = pd.read_csv(tsvfile,sep='\t')
        x = np.array(df['X'])
        xdiff = difference(np.max(x),np.min(x))
	y = np.array(df['Y'])
	ydiff = difference(np.max(y),np.min(y))
	z = np.array(df['Z'])
        zdiff = difference(np.max(z),np.min(z))	

        return xdiff, ydiff, zdiff

def difference(dim_max,dim_min):
        if float(dim_max) == float(dim_min):
                return 0
        if float(dim_max) > 0 and float(dim_min) > 0:
                return float(dim_max) - float(dim_min)
        elif float(dim_max) > 0 and float(dim_min) < 0:
                return float(dim_max) + np.abs(float(dim_min))
        elif float(dim_max) < 0 and float(dim_min) < 0:
                return float(dim_max) - float(dim_min)

proj_dir = '/oak/stanford/groups/menon/fMRI_PREP_MAIDEN_VOYAGE/fmriprep/'

# subs is a list of the subjects present in the fmriprep directory
# tasks is a list of the tasks within that subject
subs = []
tasks = set()

for dirlist, subdirs, flist in os.walk(proj_dir):
	subdir_list = subdirs
	fl_list = flist
	for element in subdir_list:
		if 'sub' in str(element):
			element = element[4:]
			subs.append(element)

	for element in fl_list:
		if 'task' in str(element) and '.nii' in str(element):
			element = element.split('-')
			t = str(element[2])
			tasks.add(t[:-11])

# Construct the subject-wise .csv
with open(proj_dir+'movement_exclusion.csv', 'wb') as csvfile:
	spamwriter = csv.writer(csvfile, delimiter=',')
	spamwriter.writerow(['task','PID','x_range','y_range','z_range','roll_range','pitch_range','yaw_range'])

	for task in list(tasks):
		for su in subs:
			try_path = proj_dir+'/sub-%s/func/sub-%s_task-%s_bold_confounds.tsv'%(su,su,task)
			print(try_path)
			if os.path.exists(try_path):
				print(try_path)
				[x_range,y_range,z_range] = extract_stats(try_path)
				print(task)
				print(su)
				print(x_range)
				print(y_range)
				print(z_range)
				spamwriter.writerow([task,su,x_range,y_range,z_range,'null','null','null'])
			else:
				continue

