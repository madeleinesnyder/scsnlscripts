
'''
    Saves a csv of which subjects to exclude based on global lab threshold
    from variance explained by motion components after running AROMA

    Fill out variables at top of script and then run:

    $ bash
    $ LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/kochalka/lib64/
    $ python AROMA_exclusion_threshold.py

'''

# Threshold of variance explained by signal to exlude by
threshold = 0.507 #Default as determined by Q1 of signal distribution on 538 task runs of lab data

# AROMA version (nonaggr or aggr)
which_aroma = 'nonaggr'

# Path to subjectlist file
subjectlist = 'test_subject_list.txt'

# List of sessions you are running for. Do not need to be present for all subjects
sessions = ['Encoding_1','encoding_1']

# Path/name of output csv for list of subjects/session to exclude
output_file = 'test_exclusion_output.csv'


import os, csv
import numpy as np
#import pandas as pd

def load_subjects(filename):
    #data = pd.read_csv(filename,header=None)
    data = np.loadtxt(filename,dtype=str)    
    subjects = data
    #subjects = data[0].values
    return subjects

subjects = load_subjects(subjectlist)
tasks = sessions

to_exclude = {'subject':[],'task':[]}
for s in range(len(subjects)):
    for t in range(len(tasks)):
        subject = subjects[s].strip("b''")
        year = '20'+subject[0:2]
        path = '/mnt/musk1/%s/%s/fmri/%s/AROMA_%s/'%(year,subject,tasks[t],which_aroma)
        #print(path)
        try:
            motion_ics = np.loadtxt(path+'classified_motion_ICs.txt',delimiter=',')
            var_exp = np.loadtxt(path+'melodic.ica/eigenvalues_percent')
            ic_indx = [m-1 for m in motion_ics]
            new_var_exp = []
            for v in range(len(var_exp)):
                if v != 0:
                    new_var_exp.append(var_exp[v] - var_exp[v-1])
                else:
                    new_var_exp.append(var_exp[v])
            var_amt = [new_var_exp[int(m)] for m in ic_indx]
            var_amt = np.sum(var_amt)

            if var_amt >= threshold:
                violations = [True]
            else:
                violations = [False]

            if True in violations:
                print('Exclude subject %s on session %s'%(subject,tasks[t]))
                to_exclude['subject'].append(subject)
                to_exclude['task'].append(tasks[t])
        except:
            pass
print('Done')

with open(output_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(to_exclude.keys())
    writer.writerows(zip(*to_exclude.values()))


#exclusion_output = pd.DataFrame.from_dict(to_exclude)
#exclusion_output.to_csv(output_file,index=False)
