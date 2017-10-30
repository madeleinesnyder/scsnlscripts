import os,sys
import csv
import pandas as pd
import numpy as np
from subprocess import call
import subprocess
import pdb
import os.path as op

script = sys.argv[0]
raw_dir = sys.argv[1]
output_dir = sys.argv[2]
names = sys.argv[3].split(',')
onsetColname = sys.argv[4]
offsetColname = sys.argv[5]
if offsetColname == 'none':
	offsetColname == ''
durcolname = sys.argv[6]
if durcolname == 'none':
	durcolname = ''
presentDur = sys.argv[7]
tasktypecolnames = sys.argv[8].split(',')
subjectlist = sys.argv[9]
runlist = sys.argv[10]
sub_df = pd.read_csv(subjectlist)
run_df = []
dummycolname = 'DummyDelay.OnsetTime'

with open(runlist) as F:
	for line in F:
		line = line[:-1]
		run_df.append(line)

# Stim_stop goSig.Acc function
def boolfunc(a,b):
        if a == 0 and b == 1:
                return 'go'
        elif a == 0 and b == 0:
                return 'trash'
        elif a == 1 and b == 0:
                return 'US'
        elif a == 1 and b == 1:
                return 'SS'

# Go through the directory and look for all the tabfiles and txt files that need to be converted
for row in sub_df.iterrows():
	for run in run_df:
		srflist,tabs = [],[]
		srdir =  raw_dir+str(row[1]['PID'])+'/visit'+str(row[1]['visit'])+'/session'+str(row[1]['session'])+'/fmri/'+run+'/behavioral/'
        	for dirlist,subdirlist,flist in os.walk(srdir):
			srflist.append(flist)
		for File in srflist:
			for element in File:
				if element[-3:] == 'txt':
					if element[0:7] == 'tabfile':
						try:
							call(['/home/groups/menon/lab_shared/software/.gem/ruby/2.4.0/bin/eprime2tabfile -o %s %s'%(srdir+element[:-3]+'tsv',srdir+element)],shell=True)
						except:
							call(['mv '+srdir+element+' '+srdir+element[:-3]+'tsv'],shell=True)
					else:
						try:
							 call(['/home/groups/menon/lab_shared/software/.gem/ruby/2.4.0/bin/eprime2tabfile -o %s %s'%(srdir+'tabfile_'+element[:-3]+'tsv',srdir+element)],shell=True)
						except:
							print("Couldn't convert file")
				elif element[-3:] == 'tsv':
                                        if element[0:7] == 'tabfile':
                                                tabs.append(element)
					else:
						call(['mv '+srdir+element+' '+srdir+'tabfile_'+element],shell=True)
		# for each tabfile, make a list of the things for the 
		for tabfile in tabs:
			# Append a new column for the tasktype
			df = pd.read_csv(srdir+tabfile,sep='\t')
			stim_stop = tasktypecolnames[0]
			resp = tasktypecolnames[1]
			df['TrialType'] = df.apply(lambda row: boolfunc(row[stim_stop],row[resp]),axis=1) 

			class nameThing:
				def __init__(self,n,onsetTimes,offsetTimes,durations):
					self.n = n
					self.onsetTimes = onsetTimes
					self.offsetTimes = offsetTimes
					self.durations = durations

			# For each run name, find all the onset times and the durations, and put them in the name class		
			taskfile = open(output_dir+'/data/imaging/participants/'+str(row[1]['PID'])+'/visit'+str(row[1]['visit'])+'/session'+str(row[1]['session'])+'/fmri/'+run+'/task_design_Test.m','w')
			print("Creating output file")
			taskfile.write('run_name = '+"'"+run+"'"+'\n')
        		taskfile.write('\n')
			# Create the nameThing for each name and append it to the taskfile
			for i,name in enumerate(names):
				obj_name = name+'_run'
				onsets,offsets,durs = [],[],[]
		
				for trial in df.iterrows():
					# Normalize the onset times associated with that trial name 
					if df[dummycolname][0] != '':
						dummydelay = df[dummycolname][0] + 9900	
						first_value = 0
					else:
						dummydelay = 0
						first_value = df[onsetColname][0] + 300
					
					# Append the onset times associtaes with that trial name
					if trial[1]['TrialType'] == name:
						onset_vec = trial[1][onsetColname] + dummydelay - first_value
						onsets.append(onset_vec/1000)
						# Append the offset time associated with that trial name (if it exists)
						try:
							offsets.append(trial[1][offsetColname])
						except:
							pass
						# Append the druation time associated with tha trial name (if it exists)
						try:
							durs.append(trial[1][durcolname])
						except:
							pass

				# Put the onset times and durations in the run object named for that trial 
				obj_name = nameThing(name,onsets,offsets,durs)
				if obj_name.durations == [] and obj_name.offsetTimes == []:
					tmp = list([float(presentDur)]*len(obj_name.onsetTimes))
					obj_name.durations = tmp
				elif obj_name.offsetTimes != []:
					tmp = [obj_name.onsetTimes[0][z]-obj_name.offsetTimes[z] for z in range(len(obj_name.onsetTimes[0]))]
					obj_name.durations = tmp
				else:
					print("Error: no offset times, present duration, or durations specified")
				# Write the file using the lists created above from the dataframe
				taskfile.write("names{"+str(i+1)+"} = '"+names[i]+"';\n")
				taskfile.write("onsets{"+str(i+1)+"} = "+str(obj_name.onsetTimes)+";\n")
				taskfile.write("durations{"+str(i+1)+"} = "+str(obj_name.durations)+";\n")
				taskfile.write("\n")
			taskfile.write('rest_exists = 1;'+"\n")
			taskfile.write("save task_design.mat sess_name names onsets durations rest_exists"+"\n")
			taskfile.close()
