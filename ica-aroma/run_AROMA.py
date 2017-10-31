import os
import os.path
import sys
import numpy as np
import scipy as sp

def source_fsl5():
  try:
    os.system('setenv FSLDIR /home/groups/menon/lab_shared/software/fsl')
    os.system('source ${FSLDIR}/etc/fslconf/fsl.csh')
    os.system('setenv PATH ${FSLDIR}/bin:${PATH}')
  except:
    print('Cannot source fsl')

def ReadSubjectList(filename):
  try:
    f = open(filename, 'rU')
    subjectIDs = f.read()
    f.close()
    subjectlist = subjectIDs.split()
  except:
    print('Cannot read subject list')
  return subjectlist

def Run_ICA_AROMA_nonaggr(input_file, output_folder, mc_file, brain_mask_fname):
  try:
    cmd = 'python /oak/stanford/groups/menon/scsnlscripts/ica-aroma/ICA-AROMA-master/ICA_AROMA.py -in ' + input_file + ' -out ' + output_folder + ' -mc ' + mc_file + ' -tr 2 ' + ' -m ' + brain_mask_fname
    print(cmd)
    os.system(cmd)
  except :
    print('Cannot run AROMA nonaggr')

def Run_ICA_AROMA_aggr(input_file, output_folder, mc_file, brain_mask_fname):
  try:
    cmd = 'python /oak/stanford/groups/menon/scsnlscripts/ica-aroma/ICA-AROMA-master/ICA_AROMA.py -in ' + input_file + ' -out ' + output_folder + ' -mc ' + mc_file + ' -tr 2 ' + ' -m ' + brain_mask_fname + ' -den aggr'
    print(cmd)
    os.system(cmd)
  except:
    print('Cannot run AROMA aggr')

def main(subject_path,brainmask,pipeline,aggr_sign):
    source_fsl5 # aggr needs fsl5
    print(subject_path)
    brain_mask_fname = brainmask
    interim_path_src = subject_path+'/smoothed_spm8/'
    interim_path_des = subject_path+'/AROMA'
    tmp1_input_file = interim_path_src+pipeline+'I.nii.gz'
    tmp2_input_file = interim_path_src+pipeline+'I.nii'

    if os.path.isfile(tmp2_input_file):
        input_file = tmp2_input_file
    elif os.path.isfile(tmp1_input_file):
        input_file = tmp1_input_file
    else:
	print(tmp1_input_file)
        print('Error: File path does not exist')
        exit()
    output_folder = interim_path_des+'_'+aggr_sign
    mc_file = interim_path_src + 'rp_I.txt'

    if aggr_sign == 'nonaggr':
        Run_ICA_AROMA_nonaggr(input_file, output_folder, mc_file, brain_mask_fname)
    elif aggr_sing == 'aggr':
        Run_ICA_AROMA_nonaggr(input_file, output_folder, mc_file, brain_mask_fname)
    else:
        print('Error: Aggr sign does not match')
        exit()



if __name__ == '__main__':
  main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

