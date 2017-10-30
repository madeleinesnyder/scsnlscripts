from subprocess import call
from scipy.io import savemat
import numpy as np
import subprocess
from numpy import array
import os
import sys
import pandas as pd
import csv
from time import sleep
import pdb

oakpath = os.environ['OAK']

class Config(object):
	'''Basic config class

	Attributes:
		subject: A subject to run your script on
		runlist: A list of runs to run your script on
		data_dir: Directory holding raw data
	'''
	data_dir = oakpath + 'rawdata'
	# self assign oakpath for bash template

	def __init__(self, subject,runlist):
		self.subject = subject
		self.runlist = runlist


	def set_data_dir(self,raw_dir):
		'''
		  Allow user to override inherited raw data directory
		'''
		self.data_dir = raw_dir

class fMRIprep(Config):
	'''
	fmri prep pipeline config class
	'''

	def __init__(self,subject,raw_dir,runlist,subjectlist,project_dir):
		super(fMRIprep, self).__init__(subject,runlist)
		self.subjectlist = subjectlist
		self.project_dir = project_dir
		self.raw_dir = raw_dir

class DTIprep(Config):
        '''
        dti prep pipeline config class
        '''

        def __init__(self,subject,raw_dir,runlist,subjectlist,project_dir):
                super(DTIprep, self).__init__(subject,runlist)
                self.subjectlist = subjectlist
                self.project_dir = project_dir
                self.raw_dir = raw_dir
		
class mriqc(Config):
        '''
        mriqc pipeline config class
        '''

        def __init__(self,subject,raw_dir,runlist,subjectlist,project_dir,working_dir):
                super(mriqc, self).__init__(subject,runlist)
                self.subjectlist = subjectlist
                self.project_dir = project_dir
	        self.working_dir = working_dir
                self.raw_dir = raw_dir

class Preprocessing(Config):
	'''Preprocessing config class

	Attributes:
	  Required:
		- subject index
        - subject list
		- runlist
		- pipeline
		- project directory
	  Optional:
		- inputimgprefix (defaults to None)
		- SPGRsubjectlist (defaults to None)
		- SPGRfilename (defaults to 'watershed_spgr')
		- tr_val (defaults to 2.0)
	'''

	def __init__(self,subject,runlist,subjectlist,pipeline,project_dir,
				 inputimgprefix='',SPGRsubjectlist='',
				 SPGRfilename='watershed_spgr',tr_val=2.0):
		super(Preprocessing, self).__init__(subject,runlist)
		self.subjectlist = subjectlist
		self.pipeline = pipeline
		self.project_dir = project_dir
		self.inputimgprefix = inputimgprefix
		self.SPGRsubjectlist = SPGRsubjectlist
		self.SPGRfilename = SPGRfilename
		self.tr_val = tr_val
		self.data_dir = self.data_dir

		# error handling

		pipeline_options = ['swar','swavr','swaor',
				    'swgcar','swgcavr','swgcaor',
				    'swfar','swfavr','swfaor',
				    'swgcfar','swgcfavr','swgcfaor']
		if self.pipeline not in pipeline_options:
			raise Exception('Incorrect Pipeline specified')

		if not os.path.isdir(self.project_dir):
			raise Exception('Project directory does not exist')


class IndividualStats(Config):
	'''IndividualStats config class


	'''

	def __init__(self,subject,subjectlist,runlist,pipeline_type,project_dir,taskdesign_file,
		     contrast_file,movement_include,fmri_type,stats_dir='stats_spm8',
		     preprocessed_dir='smoothed_spm8',tr = 2,parent_dir='',image_type='nii',
		     volrepair_include=0,volpipeline_type='',volrepaired_dir='',repairedstats_dir='',
		     template_dir = oakpath + '/scsnlscripts/spm/BatchTemplates'):
		super(IndividualStats, self).__init__(subject,runlist)
		self.subjectlist = subjectlist
		self.data_dir = self.data_dir
		self.pipeline_type = pipeline_type
		self.project_dir = project_dir
		self.taskdesign_file = taskdesign_file
		self.contrast_file = contrast_file
		self.movement_include = movement_include
		self.fmri_type = fmri_type
		self.parent_dir = parent_dir
		self.stats_dir = stats_dir
		self.preprocessed_dir = preprocessed_dir
		self.image_type = image_type
		self.tr = tr
		self.volpipeline_type = volpipeline_type
		self.volrepaired_dir = volrepaired_dir
		self.repairedstats_dir = repairedstats_dir
		self.volrepair_include = volrepair_include
		self.template_dir = template_dir
		
		# error handling

		pipeline_options = ['swar','swavr','swaor',
				    'swgcar','swgcavr','swgcaor',
				    'swfar','swfavr','swfaor',
				    'swgcfar','swgcfavr','swgcfaor']
		if self.pipeline_type not in pipeline_options:
			raise Exception('Incorrect Pipeline specified')

		if not os.path.isdir(self.project_dir):
			raise Exception('Project directory does not exist')

class MultiIndividualStats(Config):
	'''MultiIndividualStats config class

		TODO: Logic for asserting subject lists are the same length
	'''

	def __init__(self,subject,subjectlist,runlist,pipeline_type,project_dir,taskdesign_file,
		     contrast_file,movement_include,fmri_type,stats_dir='stats_spm8',
		     preprocessed_dir='smoothed_spm8',image_type='nii',
		     volrepair_include=0,volpipeline_type='',volrepaired_dir='',repairedstats_dir='',
		     template_dir = oakpath +'scsnlscripts/spm/BatchTemplates'):
		super(MultiIndividualStats, self).__init__(subject,runlist)
		self.data_dir = self.data_dir
		self.subjectlist = array(subjectlist,dtype=np.object)
		self.pipeline_type = pipeline_type
		self.project_dir = project_dir
		self.taskdesign_file = taskdesign_file
		self.contrast_file = contrast_file
		self.movement_include = movement_include
		self.fmri_type = fmri_type
		self.stats_dir = stats_dir
		self.preprocessed_dir = preprocessed_dir
		self.image_type = image_type
		self.volpipeline_type = volpipeline_type
		self.volrepaired_dir = volrepaired_dir
		self.repairedstats_dir = repairedstats_dir
		self.volrepair_include = volrepair_include
		self.template_dir = template_dir


class FunctionalConnectivity(Config):
	'''FunctionalConnectivity config class

	'''

	def __init__(self,subject,subjectlist,runlist,project_dir,pipeline,TR,
		     roi_dir,roilist,fmri_type,stats_dir,preprocessed_dir,
		     numtrunc=[0,0],data_type='nii',bandpass=0,fl=0,fh=0.2):
		super(FunctionalConnectivity, self).__init__(subject,runlist)
		self.subjectlist = subjectlist
		self.data_dir = self.data_dir
		self.project_dir = project_dir
		self.pipeline = pipeline
		self.TR = TR
		self.roi_dir = roi_dir
		self.roilist = roilist
		self.fmri_type = fmri_type
		self.stats_dir = stats_dir
		self.preprocessed_dir = preprocessed_dir
		self.numtrunc = numtrunc
		self.data_type = data_type
		self.bandpass = bandpass
		self.fl = fl
		self.fh = fh

class PPI(object):
	'''PPI config class
	'''

	def __init__(self, subject, project_dir, processed_dir, parent_dir,
		     stats_dir, roi_type, roicenter_value, roiradius_value, 
		     session_value, contrast_type, contrastweights_value):

		self.subject = subject
		self.project_dir = project_dir 
		self.processed_dir = processed_dir 
		self.parent_dir = parent_dir
		self.stats_dir = stats_dir
		self.roi_type = roi_type
		self.roicenter_value = roicenter_value
	        self.roiradius_value = roiradius_value 
		self.session_value = session_value 
		self.contrast_type = contrast_type 
		self.contrastweights_value = contrastweights_value

class PPIVolRep(object):
	'''PPIVolRep config class
	'''

	def __init__(self, project_dir, processed_dir, pipeline_type,
		     stats_dir, roi_type, roicenter_value, roiradius_value, 
		     session_value, contrast_type, contrastweights_value):

		self.subject = subject
		self.project_dir = project_dir
		self.processed_dir = processed_dir
	        self.pipeline_type = pipeline_type
	        self.stats_dir = stats_dir
	        self.roi_type = roi_type
	        self.roicenter_value = roicenter_value 
	        self.roiradius_value = roiradius_value 
	        self.session_value = session_value 
	        self.contrast_type = contrast_type
	        self.contrastweights_value = contrastweights_value

class PHIPI(object):
	'''PHIPI config class
	'''

	def __init__(self, subject, project_dir, processed_dir, parent_dir,
	             stats_dir, roiname_list, roicenter_list, roiradius_list,
	             session_value, contrast_type, contrastweights_value):
		
		self.subject = subject
		self.project_dir = project_dir
	        self.processed_dir = processed_dir
	        self.parent_dir = parent_dir
	        self.stats_dir = stats_dir
	        self.roiname_list = roiname_list
	        self.roicenter_list = roicenter_list 
	        self.roiradius_list = roiradius_list
	        self.session_value = session_value
	        self.contrast_type = contrast_type
	        self.contrastweights_value = contrastweights_value 		

class gPPI(object):
	'''gPPI config class


	'''

	def __init__(self,subject,subjectlist,project_dir,stats_dir,pipeline,
		     roifile_list,roiname_list,tasks_to_include,
		     parent_dir='',confound_names=['R1','R2','R3','R4','R5','R6']):
		self.subject = subject
		self.subjectlist = subjectlist
		self.project_dir = project_dir
		self.parent_dir = parent_dir
		self.stats_dir = stats_dir
		self.pipeline = pipeline
		self.roifile_list = roifile_list
		self.roiname_list = roiname_list
		self.tasks_to_include = tasks_to_include
		self.confound_names = confound_names


class GroupStats(object):
	'''GroupStats config class


	'''

	def __init__(self,subjectlist,output_dir,stats_dir,project_dir,
		     fmri_type,parent_dir='',regfile='',
		     template_dir = oakpath +'scsnlscripts/spm/BatchTemplates'):
		self.subjectlist = subjectlist
		self.project_dir = project_dir
		self.fmri_type = fmri_type
		self.parent_dir = parent_dir
		self.template_dir = template_dir
		self.regfile = regfile
		self.output_dir = output_dir
		self.stats_dir = stats_dir




class TaskDesign(object):
	'''TaskDesign config class

	Attributes:
	 Required:
	   - subject
	   - runlist
	   - edat_dir
	   - onsetcolumn
	   - offsetcolumn
	   - taskname
	   - matchtaskname
	   - taskdesign
	   - restexist
	'''

	def __init__(self,subject,subjectlist,runlist,edat_dir,
		     onsetcolumn,restexist,taskname,raw_dir, project_dir,
		     matchtaskname,taskdesignfilename,offsetcolumn,
		     globalonsetcolumn,globalonsetshift,
		     durationcolumn,presetduration,tasktypecolumn,runname):
		self.subject = subject
		self.subjectlist = subjectlist
		self.runlist = runlist
		self.raw_dir = raw_dir
		self.project_dir = project_dir
		self.edat_dir = edat_dir
		self.onsetcolumn = onsetcolumn
		self.offsetcolumn = offsetcolumn
		self.taskname = array(taskname, dtype=np.object)
		self.matchtaskname = array(matchtaskname, dtype=np.object)
		self.taskdesignfilename = taskdesignfilename
		self.restexist = restexist
		self.globalonsetcolumn = globalonsetcolumn
		self.globalonsetshift = globalonsetshift
		self.durationcolumn = durationcolumn
		self.presetduration = presetduration
		self.tasktypecolumn = array(tasktypecolumn, dtype=np.object)
		self.runname = runname

class RSAIndividual(object):

	def __init__(self,subject,roi_folder,roi_list,stats_folder,
	             map_type,map_pair_index,project_dir,output_dir,user_fname,marsbar_path):
		self.subject = subject
	        self.roi_folder = roi_folder
		if type(roi_list) == list:
			self.roi_list = array(roi_list)
		else:
			self.roi_list = roi_list
	        self.roi_list = roi_list
	        self.stats_folder = stats_folder
	        self.map_type = map_type
	        self.map_pair_index = array(map_pair_index)
	        self.project_dir = project_dir
	        self.output_dir = output_dir
	        self.user_fname = user_fname
	        self.marsbar_path = marsbar_path


class RSAWholeBrain(object):
	
	def __init__(self,project_dir,subject,statsfolder,maptype,mapindex,maskfile,output_dir,searchshape,searchradius):
		self.project_dir = project_dir
	        self.subject = subject
	        self.statsfolder = statsfolder
	        self.maptype = maptype
	        self.mapindex = array(mapindex)
	        self.maskfile = maskfile
	        self.output_dir = output_dir
		self.searchshape = searchshape
		self.searchradius = searchradius

class RSAGroup(object):

	def __init__(self,group_path,roi_folder,roi_list,map_folder,map_type,
	             map_index,project_dir,output_dir,user_fname,marsbar_path):

		self.group_path = group_path
	        self.roi_folder = roi_folder
		if type(roi_list) == list:
	        	self.roi_list = array(roi_list)
		else:
			self.roi_list = roi_list
	        self.map_folder = map_folder
	        self.map_type = map_type
	        self.map_index = map_index
	        self.project_dir = project_dir
	        self.output_dir = output_dir
	        self.user_fname = user_fname
	        self.marsbar_path = marsbar_path


class ContrastChange(object):
	'''ContrastChange config class
	'''

	def __init__(self, subject, project_dir,
		     contrastfile, statsdir,
		     template_dir= oakpath + 'scsnlscripts/spm/BatchTemplates'):
		self.subject = subject
		self.project_dir = project_dir
		self.contrastfile = contrastfile
		self.statsdir = statsdir
		self.template_dir = template_dir

class ContrastGenerator(object):
	'''ContrastGenerator config class
	'''

	def __init__(self, project_dir, numcontrasts, numruns, comparewithin, 
                   movementcorrection, contrastnames, contrast):
		self.project_dir = project_dir
		self.numcontrasts = numcontrasts
		self.numruns = numruns
		self.comparewithin = comparewithin
		self.movementcorrection = movementcorrection
		self.contrastnames = array(contrastnames, dtype=np.object)
		self.contrast = array(contrast, dtype=np.object)

class MovementExclusion(object):
	'''MovementExclusion config class
	'''
	
	def __init__(self, subject, runlist, raw_dir, project_dir,preprocessedfolder,
			scantoscancrit):
		self.subject = subject
		self.runlist = runlist
		self.raw_dir = raw_dir
		self.project_dir = project_dir
		self.preprocessedfolder = preprocessedfolder
		self.scantoscancrit = scantoscancrit

class ROISignalLevel(object):

	
	def __init__(self, subject, subjectlist, project_dir,parent_dir,roi_dir,roilist,stats_dir,
		     tscore_threshold,roi_result_dir,marsbar_path):
		self.subject = subject,
	        self.subjectlist = subjectlist;
                self.parent_dir = array(parent_dir),
	        self.roi_dir = roi_dir,
	        self.roilist = roilist,
	        self.stats_dir = stats_dir,
	        self.tscore_threshold = tscore_threshold,
	        self.project_dir = project_dir,
	        self.roi_result_dir = roi_result_dir,
	        self.marsbar_path = marsbar_path

class DicomToSpgr(object):
	'''DicomToSpgr config class
	'''

	def __init__(self, subject, project_dir, raw_dir, numdicom):
		self.subject = subject
		self.project_dir = project_dir
		self.raw_dir = raw_dir
		self.numdicom = numdicom

class DTIACPC(object):

	def __init__(self,subject,project_dir,raw_dir,parent_folder,participant_path,type_folder,vista_path):
		self.subject = subject
		self.project_dir = project_dir
		self.raw_dir = raw_dir
		self.parent_folder = parent_folder
		self.participant_path = participant_path
		self.type_folder = type_folder
		self.vista_path = vista_path

class DTIDCM2NII(object):

	def __init__(self,subject,project_dir,raw_dir,parent_folder,numdicom,participant_path,type_folder,vist_path):
		self.subject = subject
		self.project_dir = project_dir
		self.raw_dir = raw_dir
		self.parent_folder = parent_folder
		self.numdicom = numdicom
		self.participant_path = participant_path
		self.type_folder = type_folder
		self.vista_path = vista_path

class DTIPreprocess(object):

	def __init__(self,subject,project_dir,participant_path,outname,vista_path,outBaseName,rawprepreocess_flag,
	             rawfittensor_flag,clobber,assetflag,numBootStrapSamples,eddyCorrect,excludeVols,bsplineInterpFlag):
		self.subject = subject,
		self.project_dir = project_dir
	        self.paticipant_path = participant_path,
	        self.outname = outname,
	        self.vista_path = vista_path,
	        self.outbaseName = outBaseName,
	        self.rawpreprocess_flag = rawprepreocess_flag,
	        self.rawfittensor_flag = rawfittensor_flag,
	        self.clobber = clobber,
 	        self.assetflag = assetflag,
	        self.numBootStrapSamples = numBootStrapSamples,
	        self.eddyCorrect = eddyCorrect,
	        self.excludeVols = excludeVols,
	        self.bsplineInterpFlag = bsplineInterpFlag

class DTIBetweenROITrack(object):
	
	def __init__(self,project_dir,raw_dir,parent_folder,subject,vista_path,
	             type_folder,roilist,output_folder):
		self.project_dir = project_dir
	        self.raw_dir = raw_dir
	        self.parent_folder = parent_folder
	        self.subject = subject
		self.vista_path = vista_path
	        self.roilist = roi_list
	        self.type_folder = type_folder
	        self.output_folder = output_folder

class DTINiftiROI2DTIROI(object):

	def __init__(self,project_dir,raw_dir,parent_folder,subject,vista_path,rois,type_folder):
		self.project_dir = project_dir
	        self.raw_dir = raw_dir
	        self.parent_folder = parent_folder
	        self.subject = subject
	        self.vista_path = vista_path
	        self.rois = rois
	        self.type_folder = type_folder

class DTITrackStats(object):

	def __init__(self,subject,project_dir,raw_dir,parent_folder,subjectlist,
		     vista_path,roilist,type_folder,output_folder):
		self.subject = subject
		self.project_dir = project_dir
	        self.raw_dir = raw_dir
	        self.parent_folder = parent_folder
	        self.vista_path = vista_path
	        self.roilist = roilist
	        self.type_folder = type_folder
	        self.output_folder = output_folder

#class dti_wholebraintrack(object):

#	def __init__(self,subject,project_dir,raw_dir,parent_folder,vista_path,type_folder,
#	             preprocessed_folder,faThresh,stepSizeMm,faThreshopts,lengthThreshMm,
#	             angleThresh,wPuncture,whichAlgorithm,whichInterp,seedVoxelOffsets):


def submit_job(config_path,script,nameID,submit_script,slurm_args,config_args):
	if script == 'fmriprep':
		if config_args['freesurfer'] == 'no':
			freesurfer = '--no-freesurfer'
		elif config_args['freesurfer'] == 'yes':
			freesurfer = '' 
		if config_args['skull_strip_template'] == 'scsnl':
			sk_temp = '--skull-strip-template scsnl'
		elif config_args['skull_strip_template'] == 'poldrack_oasis':
			sk_temp = '--skull-strip-template oasis'
		if config_args['aroma'] == 'yes':
			aroma = '--use-aroma'
		elif config_args['aroma'] == 'no':
			aroma = ''
		if config_args['smoothing'] == 4:
			smoothing = '--smoothing-kernel 4'
		elif config_args['smoothing'] == 6:
			smoothing = '--smoothing-kernel 6'
                elif config_args['smoothing'] == 8:
			smoothing = '--smoothing-kernel 8'
		elif config_args['smoothing'] == 10:
                        smoothing = '--smoothing-kernel 10'
		if config_args['pipeline'] == 'swar':
			coreg_temp = '--output-space template'
		elif config_args['pipeline'] == 'swgcar':
			coreg_temp = '--output-space T1w'
		if config_args['task_t'] == 'all':
			task_t = ''
		else:
			task_t = '-t '+config_args['task_t']
		if config_args['anat_only'] == 'no':
			anat_only = ''
		elif config_args['anat_only'] == 'yes':
			anat_only = '--anat-only'
                if config_args['func_only'] == 'yes':
			func_only = '--dismiss-t1w'
			submit_script = '/oak/stanford/groups/menon/scsnlscripts/spm/Utility/submit_singularityfmriprep_experimental.sh'
		elif config_args['func_only'] == 'no':
			func_only = ''
		nameID = nameID.split('_')
		Pnum = nameID[0]
		Vnum = nameID[1][-1:]
		Snum = nameID[2][-1:]
		nameID = nameID[0]+' '+Vnum+' '+Snum
		callstr = ('ml biology fsl && sbatch'
			   ' -J %s'
			   ' -o %s'
			   ' -e %s'
			   ' -t %s'
			   ' --qos=%s'
			   ' -N %s'
			   ' --mem-per-cpu=%s'
			   ' -p %s'
			   ' %s "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s"')%(slurm_args['job_name'],
								      slurm_args['job_output'],
								      slurm_args['job_error'],
								      slurm_args['job_time'],
								      slurm_args['job_qos'],
								      slurm_args['job_nnodes'],
								      slurm_args['job_mem'],
								      slurm_args['job_partition'],
								      submit_script,
							 	      config_args['raw_dir'],
								      config_args['project_dir'],
								      nameID,
								      Pnum,
								      Vnum,
								      Snum,
								      config_args['project_dir']+config_args['working_dir'],
								      aroma,
								      freesurfer,
								      sk_temp,
								      coreg_temp,
								      task_t,
								      smoothing,
								      anat_only,
 								      func_only)
		call(callstr,shell=True)

	elif script == 'dtiprep':	
		nameID = nameID.split('_')
                Pnum = nameID[0]
                Vnum = nameID[1][-1:]
                Snum = nameID[2][-1:]
		callstr = ('sbatch'
			   ' -J %s'
			   ' -o %s'
			   ' -e %s'
			   ' -t %s'
			   ' --qos=%s'
			   ' -N %s'
			   ' --mem-per-cpu=%s'
			   ' -p %s'
			   ' %s "%s" "%s" "%s" "%s"')%(slurm_args['job_name'],
                                                                      slurm_args['job_output'],
                                                                      slurm_args['job_error'],
                                                                      slurm_args['job_time'],
                                                                      slurm_args['job_qos'],
                                                                      slurm_args['job_nnodes'],
                                                                      slurm_args['job_mem'],
                                                                      slurm_args['job_partition'],
                                                                      submit_script,
                                                                      config_args['project_dir']+'data/imaging/BIDS/',
                                                                      Pnum,
                                                                      Vnum,
								      Snum)
		call(callstr,shell=True)

	elif script == 'mriqc':
                callstr = ('sbatch'
                           ' -J %s'
                           ' -o %s'
                           ' -e %s'
                           ' -t %s'
                           ' --qos=%s'
                           ' -N %s'
                           ' --mem-per-cpu=%s'
                           ' -p %s'
                           ' %s "%s" "%s" "%s"')%(slurm_args['job_name'],
                                                                      slurm_args['job_output'],
                                                                      slurm_args['job_error'],
                                                                      slurm_args['job_time'],
                                                                      slurm_args['job_qos'],
                                                                      slurm_args['job_nnodes'],
                                                                      slurm_args['job_mem'],
                                                                      slurm_args['job_partition'],
                                                                      submit_script,
                                                                      config_args['project_dir']+'data/imaging/BIDS/',
                                                                      config_args['project_dir']+'mriqc/',
                                                                      config_args['project_dir'])
                call(callstr,shell=True)

		
  	elif script == 'taskdesign_generation':
		callstr = ('sbatch'
                           ' -J %s'
                           ' -o %s'
                           ' -e %s'
                           ' -t %s'
                           ' --qos=%s'
                           ' -N %s'
                           ' --mem-per-cpu=%s'
                           ' -p %s'
                           ' %s "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s" "%s"')%(slurm_args['job_name'],
                                                  slurm_args['job_output'],
                                                  slurm_args['job_error'],
                                                  slurm_args['job_time'],
                                                  slurm_args['job_qos'],
                                                  slurm_args['job_nnodes'],
                                                  slurm_args['job_mem'],
                                                  slurm_args['job_partition'],
                                                  submit_script,
                                                  config_args['raw_dir'],
                                                  config_args['project_dir'],
						  config_args['taskname'],
						  config_args['onsetcolumn'],
						  config_args['offsetcolumn'],
						  config_args['durationcolumn'],
						  config_args['presetduration'],
						  config_args['tasktypecolumn'],
						  config_args['subjectlist'],
						  config_args['runlist'])
                call(callstr,shell=True)

	else:
                callstr = ('sbatch'
                           ' -J %s'
                           ' -o %s'
                           ' -e %s'
                           ' -t %s'
                           ' --qos=%s'
                           ' -N %s'
                           ' --mem-per-cpu=%s'
                           ' -p %s'
                           ' %s "%s(\'%s/config_args_%s_%s.mat\')"')%(slurm_args['job_name'],
                                                                      slurm_args['job_output'],
                                                                      slurm_args['job_error'],
                                                                      slurm_args['job_time'],
                                                                      slurm_args['job_qos'],
                                                                      slurm_args['job_nnodes'],
                                                                      slurm_args['job_mem'],
                                                                      slurm_args['job_partition'],
                                                                      submit_script,
                                                                      script,
                                                                      config_path,
                                                                      script,
                                                                      nameID)
                call(callstr,shell=True)


def run_script(config_args,slurm_args,script):

	scratch_dir = os.environ['SCRATCH']
	config_path = '%s/config_files'%(scratch_dir)
	submit_script = '/oak/stanford/groups/menon/scsnlscripts/spm/Utility/submit_job.sh'

	if not os.path.isdir(config_path):
		call('mkdir %s'%(config_path),shell=True)

	if script == 'rsa_group':
		config = RSAGroup(group_path = config_args['group_path'],
				  roi_folder = config_args['roi_folder'],
				  roi_list = config_args['roi_list'],
				  map_folder = config_args['map_folder'],
				  map_type = config_args['map_type'],
	             		  map_index = config_args['map_index'],
				  project_dir = config_args['project_dir'],
				  output_dir = config_args['output_dir'],
				  user_fname = config_args['user_fname'],
				  marsbar_path = config_args['marsbar_path'])
		nameID = ''
		input_args=config.__dict__
		savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
		submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)


	if script == 'contrast_generator':
		config = ContrastGenerator(project_dir = config_args['project_dir'],
					   numcontrasts = config_args['numcontrasts'],
					   numruns = config_args['numruns'],
					   comparewithin = config_args['comparewithin'],
					   movementcorrection = config_args['movementcorrection'],
					   contrastnames = config_args['contrastnames'],
					   contrast = config_args['contrast'])
		nameID = ''
		input_args=config.__dict__
		savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
		submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)

	elif 'groupstats' in script:
		subjectlist = config_args['subjectlist']

		config = GroupStats(subjectlist=subjectlist,
				    stats_dir=config_args['stats_dir'],
				    project_dir=config_args['project_dir'],
				    fmri_type=config_args['fmri_type'],
				    parent_dir=config_args['parent_dir'],
				    output_dir=config_args['output_dir'],
				    regfile=config_args['regfile'],
				    template_dir=config_args['template_dir'])


		nameID = config_args['stats_dir']
		input_args = config.__dict__
		savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
		submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)

	else:
		subjectlist = config_args['subjectlist']


		if script == 'multi_individualstats':
			nlists = len(subjectlist)
			nameID = ''
			Subjectlist = []
			for s_i in range(nlists):
				Subjectlist.append(pd.read_csv(subjectlist[s_i]))
				for subject_i, subject in Subjectlist[s_i].iterrows():
					subjectPID = str(subject['PID']).zfill(4)
					nameID+='%s_'%(subjectPID)
			nameID = nameID+'multiindividualstats'


			for subject_i, subject in pd.read_csv(subjectlist[0]).iterrows():
				subject_i = subject_i+1
				config = MultiIndividualStats(subject=subject_i,
							 subjectlist= config_args['subjectlist'],
							 image_type=config_args['image_type'],
							 pipeline_type=config_args['pipeline_type'],
							 project_dir=config_args['project_dir'],
							 runlist=config_args['runlist'],
							 taskdesign_file=config_args['taskdesign_file'],
							 contrast_file=config_args['contrast_file'],
							 preprocessed_dir=config_args['preprocessed_dir'],
							 stats_dir=config_args['stats_dir'],
							 movement_include=config_args['movement_include'],
							 volrepair_include=config_args['volrepair_include'],
							 volpipeline_type=config_args['volpipeline_type'],
							 volrepaired_dir=config_args['volrepaired_dir'],
							 repairedstats_dir=config_args['repairedstats_dir'],
							 template_dir=config_args['template_dir'],
							 fmri_type=config_args['fmri_type'])
				config.set_data_dir(config_args['raw_dir'])

				input_args = config.__dict__
				savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
				submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)
			
		else:
			subjectlist = pd.read_csv(subjectlist)
    
			for subject_i, subject in subjectlist.iterrows():
				subject_i = subject_i+1
				subjectPID = str(subject['PID']).zfill(4)
				nameID = '%s_visit%s_session%s'%(subjectPID,subject['visit'],subject['session'])
				

				if script == 'preprocessfmri':
					config = Preprocessing(subject=subject_i,
							       subjectlist=config_args['subjectlist'],
							       runlist=config_args['runlist'],
							       tr_val=config_args['tr_val'],
							       inputimgprefix=config_args['inputimgprefix'],
							       pipeline=config_args['pipeline'],
							       SPGRfilename=config_args['SPGRfilename'],
							       project_dir=config_args['project_dir'])
					config.set_data_dir(config_args['raw_dir'])

				elif script == 'individualstats':

					config = IndividualStats(subject=subject_i,
								 subjectlist=config_args['subjectlist'],
								 image_type=config_args['image_type'],
								 pipeline_type=config_args['pipeline_type'],
								 project_dir=config_args['project_dir'],
								 runlist=config_args['runlist'],
								 taskdesign_file=config_args['taskdesign_file'],
								 contrast_file=config_args['contrast_file'],
								 preprocessed_dir=config_args['preprocessed_dir'],
								 stats_dir=config_args['stats_dir'],
								 movement_include=config_args['movement_include'],
								 volrepair_include=config_args['volrepair_include'],
								 volpipeline_type=config_args['volpipeline_type'],
								 volrepaired_dir=config_args['volrepaired_dir'],
								 repairedstats_dir=config_args['repairedstats_dir'],
								 template_dir=config_args['template_dir'],
								 fmri_type=config_args['fmri_type'],
								 tr=config_args['tr'])
					config.set_data_dir(config_args['raw_dir'])

				elif script == 'movement_exclusion':

					config = MovementExclusion(subject=subject_i,
								    runlist=config_args['runlist'],
								    raw_dir=config_args['raw_dir'],
								    project_dir=config_args['project_dir'],
								    preprocessedfolder=config_args['preprocessedfolder'],
								    scantoscancrit=config_args['scantoscancrit'])

				elif script == 'fconnect':

					config = FunctionalConnectivity(subject=subject_i,
									subjectlist=config_args['subjectlist'],
									data_type=config_args['data_type'],
									project_dir=config_args['project_dir'],
									runlist=config_args['run'],
									stats_dir=config_args['stats_dir'],
									preprocessed_dir=config_args['preprocessed_dir'],
									pipeline=config_args['pipeline'],
									TR=config_args['TR'],
									bandpass=config_args['bandpass'],
									fl=config_args['fl'],
									fh=config_args['fh'],
									roi_dir=config_args['roi_dir'],
									roilist=config_args['roilist'],
									numtrunc=config_args['numtrunc'],
									fmri_type=config_args['fmri_type'])
					config.set_data_dir(config_args['raw_dir'])

				elif script == 'effconn_ppi':

					config = PPI(subject=subject_i,
						     project_dir=config_args['project_dir'], 
						     processed_dir=config_args['processed_dir'], 
						     parent_dir=config_args['parent_dir'],
						     stats_dir=config_args['stats_dir'],
						     roi_type=config_args['roi_type'],
						     roicenter_value=config_args['roicenter_value'],
						     roiradius_value=config_args['roiradius_value'],
						     session_value=config_args['session_value'], 
						     contrast_type=config_args['contrast_type'], 
						     contrastweights_value=config_args['contrastweights_value'])

				elif script == 'effconn_ppi_volrep':

					config = PPIVolRep(subject=subject_i,
							   project_dir=config_args['project_dir'],
							   processed_dir=config_args['processed_dir'],
							   pipeline_type=config_args['pipeline_type'],
							   stats_dir=config_args['stats_dir'],
							   roi_type=config_args['roi_type'],
							   roicenter_value=config_args['roicenter_value'], 
							   roiradius_value=config_args['roiradius_value'],
							   session_value=config_args['session_value'], 
							   contrast_type=config_args['contrast_type'],
							   contrastweights_value=config_args['contrastweights_value'])

				elif script == 'effconn_phipi':

					config = PHIPI(subject=subject_i,
						       project_dir=config_args['project_dir'],
						       processed_dir=config_args['processed_dir'],
						       parent_dir=config_args['parent_dir'],
						       stats_dir=config_args['stats_dir'],
						       roiname_list=config_args['roiname_list'],
						       roicenter_list=config_args['roicenter_list'], 
						       roiradius_list=config_args['roiradius_list'],
						       session_value=config_args['session_value'],
						       contrast_type=config_args['contrast_type'],
						       contrastweights_value=config_args['contrastweights_value'])

				elif script == 'scsnl_gPPI':

					config = gPPI(subject=subject_i,
						      subjectlist = config_args['subjectlist'],
						      project_dir=config_args['project_dir'],
						      parent_dir=config_args['parent_dir'],
						      stats_dir=config_args['stats_dir'],
						      pipeline=config_args['pipeline'],
						      roifile_list=config_args['roifile_list'],
						      roiname_list=config_args['roiname_list'],
						      tasks_to_include=config_args['tasks_to_include'],
						      confound_names=config_args['confound_names'])

				elif script == 'rsa_individual':

					config = RSAIndividual(subject=subject_i,
							       roi_folder=config_args['roi_folder'],
							       roi_list=config_args['roi_list'],
							       stats_folder=config_args['stats_folder'],
			     				       map_type=config_args['map_type'],
							       map_pair_index=config_args['map_pair_idex'],
							       project_dir=config_args['projectsdir'],
							       output_dir=config_args['output_dir'],
							       user_fname=config_args['user_fname'],
							       marsbar_path=config_args['marsbar_path'])


				elif script == 'rsa_wholebrain':

					config = RSAWholeBrain(subject=subject_i,
							       project_dir=config_args['project_dir'],
							       statsfolder=config_args['statsfolder'],
							       maptype=config_args['maptype'],
							       mapindex=config_args['mapindex'],
							       maskfile=config_args['maskfile'],
							       output_dir=config_args['output_dir'],
							       searchshape=config_args['searchshape'],
							       searchradius=config_args['searchradius'])

				elif script == 'roi_signallevel':

					config = ROISignalLevel(subject=subject_i,
								subjectlist=config_args['subjectlist'],
								parent_dir=config_args['parent_dir'],
								roi_dir=config_args['roi_dir'],
								roilist=config_args['roilist'],
								stats_dir=config_args['stats_dir'],
								tscore_threshold=config_args['tscore_threshold'],
								project_dir=config_args['project_dir'],
								roi_result_dir=config_args['roi_result_dir'],
								marsbar_path=config_args['marsbar_path'])

				elif script == 'contrastchange':
	
					config = ContrastChange(subject=subject_i,
								project_dir=config_args['project_dir'],
								contrastfile=config_args['contrastfile'],
								statsdir=config_args['statsdir'],
								template_dir=config_args['template_dir'])

				elif script == 'dicomtospgr':

					config = DicomToSpgr(subject=subject_i,
							     project_dir=config_args['project_dir'],
							     raw_dir=config_args['raw_dir'],
							     numdicom=config_args['numdicom'])

				elif script == 'dti_acpc':

					config = DTIACPC(subject=subject_i,
							 project_dir=config_args['project_dir'],
							 raw_dir=config_args['raw_dir'],
							 parent_folder=config_args['parent_folder'],
							 participant_path=config_args['participant_path'],
							 type_folder=config_args['type_folder'],
							 vista_path=config_args['vista_path'])

				elif script == 'dti_dcm2nii':

					config = DTIDCM2NII(subject=subject_i,
							    project_dir=config_args['project_dir'],
							    raw_dir=config_args['raw_dir'],
							    parent_folder=config_args['parent_folder'],
							    numdicom=config_args['numdicom'],
							    participant_path=config_args['participant_path'],
							    type_folder=config_args['type_folder'],
							    vista_path=config_args['vista_path'])

				elif script == 'dti_preprocess':

					config = DTIPreprocess(subject=subject_i,
							       project_dir=config_args['project_dir'],
		       					       paticipant_path= ['participant_path'],
		       					       outname=config_args['outname'],
		       					       vista_path=config_args['vista_path'],
		       					       outbaseName=config_args['outBaseName'],
		       					       rawpreprocess_flag=config_args['rawprepreocess_flag'],
		       					       rawfittensor_flag=config_args['rawfittensor_flag'],
		       					       clobber=config_args['clobber'],
	 	       					       assetflag=config_args['assetflag'],
		       					       numBootStrapSamples=config_args['numBootStrapSamples'],
		       					       eddyCorrect=config_args['eddyCorrect'],
		       					       excludeVols=['excludeVols'],
		       					       bsplineInterpFlag=['bsplineInterpFlag'])

				elif script == 'dti_betweenroitrack':

					config = DTIBetweenROITrack(subject=subject_i,
								    project_dir=config_args['project_dir'],
		       						    raw_dir=config_args['raw_dir'],
		       						    parent_folder=config_args['parent_folder'],
		       						    vista_path=config_args['vista_path'],
		       						    roilist=config_args['type_folder'],
		       						    type_folder=config_args['type_folder'],
		       						    output_folder=config_args['output_folder'])

				elif script == 'dti_niftiroi2dtiroi':

					config = DTINiftiROI2DTIROI(subject=subject_i,
							            project_dir=config_args['project_dir'],
		       						    raw_dir=config_args['raw_dir'],
		       						    parent_folder=config_args['parent_folder'],
		       						    vista_path=config_args['vista_path'],
		       						    rois=config_args['rois'],
		       						    type_folder=config_args['type_folder'])

				elif script == 'dti_trackstat':

					config = DTITrackStats(subject=subject_i,
							       project_dir=config_args['project_dir'],
		       					       raw_dir=config_args['raw_dir'],
		       					       parent_folder=config_args['parent_folder'],
		       					       vista_path=config_args['vista_path'],
		       					       roilist=config_args['roilist'],
		       					       type_folder=config_args['type_folder'],
		       					       output_folder=config_args['output_folder'])

                                elif script == 'taskdesign_generation':
					submit_script = '/oak/stanford/groups/menon/scsnlscripts/spm/Utility/submit_pyjob.sh'
	
                                        config = TaskDesign(subject=subject_i,
                                                            subjectlist=config_args['subjectlist'],
                                                            runlist=config_args['runlist'],
                                                            raw_dir=config_args['raw_dir'],
                                                            project_dir=config_args['project_dir'],
                                                            edat_dir=config_args['edat_dir'],
                                                            globalonsetcolumn=config_args['globalonsetcolumn'],
                                                            globalonsetshift=config_args['globalonsetshift'],
                                                            onsetcolumn=config_args['onsetcolumn'],
                                                            offsetcolumn=config_args['offsetcolumn'],
                                                            durationcolumn=config_args['durationcolumn'],
                                                            presetduration=config_args['presetduration'],
                                                            tasktypecolumn=config_args['tasktypecolumn'],
                                                            taskname=config_args['taskname'],
                                                            matchtaskname=config_args['matchtaskname'],
                                                            runname=config_args['runname'],
                                                            taskdesignfilename=config_args['taskdesignfilename'],
                                                            restexist=config_args['restexist'])

                                        input_args = config.__dict__

                                        submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)
                                        sleep(1)


				elif script == 'fmriprep':
				        submit_script = '/oak/stanford/groups/menon/scsnlscripts/spm/Utility/submit_singularityfmriprep.sh'
					
					config = fMRIprep(subject=subject_i,
							  raw_dir=config_args['raw_dir'],
							  subjectlist=config_args['subjectlist'],
							  runlist=config_args['runlist'],
							  project_dir=config_args['project_dir'])
					config.set_data_dir(config_args['raw_dir'])

					input_args = config.__dict__

					submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)
					sleep(1)

				elif script == 'dtiprep':
                                        submit_script = '/oak/stanford/groups/menon/scsnlscripts/spm/Utility/submit_dtijob.sh'

                                        config = DTIprep(subject=subject_i,
                                                          raw_dir=config_args['raw_dir'],
                                                          subjectlist=config_args['subjectlist'],
                                                          runlist=config_args['runlist'],
                                                          project_dir=config_args['project_dir'])
                                        config.set_data_dir(config_args['raw_dir'])

                                        input_args = config.__dict__

                                        submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)
                                        sleep(1)

			        elif script == 'mriqc':
                                        submit_script = '/oak/stanford/groups/menon/scsnlscripts/spm/Utility/submit_singularitymriqc.sh'

                                        config = mriqc(subject=subject_i,
                                                          raw_dir=config_args['raw_dir'],
                                                          subjectlist=config_args['subjectlist'],
                                                          runlist=config_args['runlist'],
                                                          project_dir=config_args['project_dir'],
							  working_dir=config_args['working_dir'])
                                        config.set_data_dir(config_args['raw_dir'])

                                        input_args = config.__dict__

                                        submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)
                                        sleep(1)

				input_args = config.__dict__
				savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
				submit_job(config_path,script,nameID,submit_script,slurm_args,config_args)
				sleep(1)
