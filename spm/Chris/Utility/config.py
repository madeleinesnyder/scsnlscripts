from subprocess import call
from scipy.io import savemat
import numpy as np
from numpy import array
import os
import sys

class Config(object):
	'''Basic config class

	Attributes:
		subject: A subject to run your script on
		sessionlist: A list of sessions to run your script on
		datadir: Directory holding raw data
	'''

	datadir='/home/groups/menon/oak/rawdata/scsnl'


	def __init__(self, subject,sessionlist):
		self.subject = subject
		self.sessionlist = sessionlist


	def set_datadir(self,rawdir):
		'''
		  Allow user to override inherited raw data directory
		'''
		self.datadir = rawdir


class Preprocessing(Config):
	'''Preprocessing config class

	Attributes:
	  Required:
		- subject
		- sessionlist
		- pipeline
		- project directory
	  Optional:
		- inputimgprefix (defaults to None)
		- SPGRsubjectlist (defaults to None)
		- SPGRfilename (defaults to 'watershed_spgr')
		- tr_val (defaults to 2.0)
	'''

	def __init__(self,subject,sessionlist,pipeline,projectdir,
				 inputimgprefix='',SPGRsubjectlist='',
				 SPGRfilename='watershed_spgr',tr_val=2.0):
		super(Preprocessing, self).__init__(subject,sessionlist)
		self.pipeline = pipeline
		self.projectdir = projectdir
		self.inputimgprefix = inputimgprefix
		self.SPGRsubjectlist = SPGRsubjectlist
		self.SPGRfilename = SPGRfilename
		self.tr_val = tr_val
		self.datadir = self.datadir

		# error handling

		pipeline_options = ['swar','swavr','swaor',
				    'swgcar','swgcavr','swgcaor',
				    'swfar','swfavr','swfaor',
				    'swgcfar','swgcfavr','swgcfaor']
		if self.pipeline not in pipeline_options:
			raise Exception('Incorrect Pipeline specified')

		if not os.path.isdir(self.projectdir):
			raise Exception('Project directory does not exist')


class IndividualStats(Config):
	'''IndividualStats config class


	'''

	def __init__(self,subject,sessionlist,pipeline_type,project_dir,taskdesign_file,
		     contrast_file,movement_include,fmri_type,stats_dir='stats_spm8',
		     preprocessed_dir='smoothed_spm8',parent_dir='',image_type='nii',
		     volrepair_include=0,volpipeline_type='',volrepaired_dir='',repairedstats_dir='',
		     template_dir='/home/groups/menon/oak/scsnlscripts/spm/BatchTemplates'):
		super(IndividualStats, self).__init__(subject,sessionlist)
		self.datadir = self.datadir
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
		if self.pipeline not in pipeline_options:
			raise Exception('Incorrect Pipeline specified')

		if not os.path.isdir(self.projectdir):
			raise Exception('Project directory does not exist')

class MultiIndividualStats(Config):
	'''MultiIndividualStats config class

		TODO: Logic for asserting subject lists are the same length
	'''

	def __init__(self,subject,sessionlist,pipeline_type,project_dir,taskdesign_file,
		     contrast_file,movement_include,fmri_type,stats_dir='stats_spm8',
		     preprocessed_dir='smoothed_spm8',image_type='nii',
		     volrepair_include=0,volpipeline_type='',volrepaired_dir='',repairedstats_dir='',
		     template_dir='/home/groups/menon/oak/scsnlscripts/spm/BatchTemplates'):
		super(MultiIndividualStats, self).__init__(subject,sessionlist)
		self.datadir = self.datadir
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

	def __init__(self,subject,sessionlist,project_dir,pipeline,TR,
		     roi_dir,roi_list,fmri_type,stats_dir,preprocessed_dir,
		     numtrunc=[0,0],data_type='nii',nonyear_dir='',bandpass=0,fl=0,fh=0.2):
		super(FunctionalConnectivity, self).__init__(subject,sessionlist)
		self.datadir = self.datadir
		self.project_dir = project_dir
		self.pipeline = pipeline
		self.TR = TR
		self.roi_dir = roi_dir
		self.roi_list = roi_list
		self.fmri_type = fmri_type
		self.stats_dir = stats_dir
		self.preprocessed_dir = preprocessed_dir
		self.numtrunc = numtrunc
		self.data_type = data_type
		self.nonyear_dir = nonyear_dir
		self.bandpass = bandpass
		self.fl = fl
		self.fh = fh

class PPI(object):
	'''PPI config class
	'''

	def __init__(self, subject, projectdir, processed_dir, parent_dir,
		     stats_dir, roi_type, roicenter_value, roiradius_value, 
		     session_value, contrast_type, contrastweights_value)

		self.subject = subject
		self.projectdir = projectdir 
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

	def __init__(self, projectdir, processed_dir, pipeline_type,
		     stats_dir, roi_type, roicenter_value, roiradius_value, 
		     session_value, contrast_type, contrastweights_value)

		self.subject = subject
		self.projectdir = projectdir
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

	def __init__(self, subject, projectdir, processed_dir, parent_dir,
	             stats_dir, roiname_list, roicenter_list, roiradius_list,
	             session_value, contrast_type, contrastweights_value)
		
		self.subject = subject
		self.projectdir = projectdir
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

	def __init__(self,subject,project_dir,stats_dir,pipeline,
		     roifile_list,roiname_list,tasks_to_include,
		     parent_dir='',confound_names=['R1','R2','R3','R4','R5','R6']):
		self.subject = subject
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

	def __init__(self,subjectlist,outputfolder,statsfolder,project_dir,
		     fmri_type,parent_dir='',regfile=[''],
		     template_dir='/home/fmri/fmrihome/SPM/spm8_scripts/BatchTemplates'):
		self.subjectlist = subjectlist
		self.project_dir = project_dir
		self.fmri_type = fmri_type
		self.parent_dir = parent_dir
		self.template_dir = template_dir
		self.regfile = regfile
		self.outputfolder = outputfolder
		self.statsfolder = statsfolder



class TaskDesign(object):
	'''TaskDesign config class

	Attributes:
	 Required:
	   - subject
	   - sessionlist
	   - edatfolder
	   - onsetcolumn
	   - offsetcolumn
	   - taskname
	   - matchtaskname
	   - taskdesign
	   - restexist
	'''

	def __init__(self,subject,sessionlist,edatfolder,
				 onsetcolumn,restexist,taskname,rawdir, projectdir,
				 matchtaskname,taskdesignfilename,offsetcolumn,
				 globalonsetcolumn,globalonsetshift,
				 durationcolumn,presetduration,
				 tasktypecolumn,sessionname):
		self.subject = subject
		self.sessionlist = sessionlist
		self.rawdir = rawdir
		self.projectdir = projectdir
		self.edatfolder = edatfolder
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
		self.sessionname = sessionname


class ContrastChange(object):
	'''ContrastChange config class
	'''

	def __init__(self, subject, projectdir,
		     contrastfile, statsdir,
		     template_dir='/home/groups/menon/oak/scsnlscripts/spm/BatchTemplates'):
		self.subject = subject
		self.projectdir = projectdir
		self.contrastfile = contrastfile
		self.statsdir = statsdir
		self.template_dir = template_dir

class ContrastGenerator(object):
	'''ContrastGenerator config class
	'''

	def __init__(self, projectdir, numcontrasts, numsessions, comparewithin, 
                   movementcorrection, contrastnames, contrast):
		self.projectdir = projectdir
		self.numcontrasts = numcontrasts
		self.numsessions = numsessions
		self.comparewithin = comparewithin
		self.movementcorrection = movementcorrection
		self.contrastnames = array(contrastnames, dtype=np.object)
		self.contrast = array(contrast, dtype=np.object)

class MovementExclusion(object):
	'''MovementExclusion config class
	'''
	
	def __init__(self,subject,sessionlist, rawdir, projectdir,preprocessedfolder,
			scantoscancrit):
		self.subject = subject
		self.sessionlist = sessionlist
		self.rawdir = rawdir
		self.projectdir = projectdir
		self.preprocessedfolder = preprocessedfolder
		self.scantoscancrit = scantoscancrit



def submit_job(config_path,script,nameID,submit_script,slurm_args):
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
	submit_script = '/home/groups/menon/oak/scsnlscripts/spm/Utility/submit_job.sh'

	call('mkdir %s'%(config_path),shell=True)

	if script == 'contrast_generator':
		config = ContrastGenerator(projectdir = config_args['projectdir'],
					   numcontrasts = config_args['numcontrasts'],
					   numsessions = config_args['numsessions'],
					   comparewithin = config_args['comparewithin'],
					   movementcorrection = config_args['movementcorrection'],
					   contrastnames = config_args['contrastnames'],
					   contrast = config_args['contrast'])
		nameID = ''
		input_args=config.__dict__
		savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
		submit_job(config_path,script,nameID,submit_script,slurm_args)


	elif 'groupstats' in script:
		subjectlist = config_args['subjectlist']
		config = GroupStats(subjectlist=subjectlist,
				    project_dir=config_args['project_dir'],
				    fmri_type=config_args['fmri_type'],
				    parent_dir=config_args['parent_dir'],
				    statsfolder=config_args['statsfolder'],
				    outputfolder=config_args['outputfolder'],
				    regfile=config_args['regfile'],
				    template_dir=config_args['template_dir'])
		nameID = config_args['statsfolder']
		input_args = config.__dict__
		savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
		submit_job(config_path,script,nameID,submit_script,slurm_args)
	else:
		subjectlist = config_args['subjectlist']
		if script == 'multi_individualstats':

			sublists = [open(sublist,"r").readlines() for sublist in subjectlist]
			n_sublists = len(sublists)

			subjectlist = []
			for isub in range(len(sublists[0])):
				sublist = [sublists[nlist][isub].replace("\n","") for nlist in range(len(sublists))]
				subjectlist.append(sublist)

		else:
			subjectfile = open(subjectlist,"r")
			subjectlist = subjectfile.readlines()

		for subject in subjectlist:
			if script != 'multi_individualstats':
				subject = subject.replace("\n","")
				nameID = subject
			else:
				nameID = '%s_%s'%(subject[0],subject[1])

			if script == 'preprocessfmri':
				config = Preprocessing(subject=subject,
						       sessionlist=config_args['sessionlist'],
						       tr_val=config_args['tr_val'],
						       inputimgprefix=config_args['inputimgprefix'],
						       pipeline=config_args['pipeline'],
						       SPGRfilename=config_args['SPGRfilename'],
						       projectdir=config_args['projectdir'])
				config.set_datadir(config_args['rawdir'])

			elif script == 'movement_exclusion':

				config = MovementExclusion(subject=subject,
							    sessionlist=config_args['sessionlist'],
							    rawdir=config_args['rawdir'],
							    projectdir=config_args['projectdir'],
							    preprocessedfolder=config_args['preprocessedfolder'],
							    scantoscancrit=config_args['scantoscancrit'])

			elif script == 'taskdesign_generation':

				config = TaskDesign(subject=subject,
						    sessionlist=config_args['sessionlist'],
						    rawdir=config_args['rawdir'],
						    projectdir=config_args['projectdir'],
						    edatfolder=config_args['edatfolder'],
						    globalonsetcolumn=config_args['globalonsetcolumn'],
						    globalonsetshift=config_args['globalonsetshift'],
						    onsetcolumn=config_args['onsetcolumn'],
						    offsetcolumn=config_args['offsetcolumn'],
						    durationcolumn=config_args['durationcolumn'],
						    presetduration=config_args['presetduration'],
						    tasktypecolumn=config_args['tasktypecolumn'],
						    taskname=config_args['taskname'],
						    matchtaskname=config_args['matchtaskname'],
						    sessionname=config_args['sessionname'],
						    taskdesignfilename=config_args['taskdesignfilename'],
						    restexist=config_args['restexist'])

			elif script == 'individualstats':

				config = IndividualStats(subject=subject,
							 image_type=config_args['image_type'],
							 pipeline_type=config_args['pipeline_type'],
							 project_dir=config_args['project_dir'],
							 parent_dir=config_args['parent_dir'],
							 sessionlist=config_args['sessionlist'],
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
				config.set_datadir(config_args['rawdir'])

			elif script == 'multi_individualstats':

				config = MultiIndividualStats(subject=subject,
							 image_type=config_args['image_type'],
							 pipeline_type=config_args['pipeline_type'],
							 project_dir=config_args['project_dir'],
							 sessionlist=config_args['sessionlist'],
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
				config.set_datadir(config_args['rawdir'])

			elif script == 'fconnect':

				config = FunctionalConnectivity(subject=subject,
							        data_type=config_args['data_type'],
							        nonyear_dir=config_args['nonyear_dir'],
							        project_dir=config_args['project_dir'],
							        sessionlist=config_args['session'],
								stats_dir=config_args['stats_dir'],
								preprocessed_dir=config_args['preprocessed_dir'],
							        pipeline=config_args['pipeline'],
							        TR=config_args['TR'],
							        bandpass=config_args['bandpass'],
							        fl=config_args['fl'],
							        fh=config_args['fh'],
							        roi_dir=config_args['roi_dir'],
							        roi_list=config_args['roi_list'],
							        numtrunc=config_args['numtrunc'],
							        fmri_type=config_args['fmri_type'])
				config.set_datadir(config_args['rawdir'])

			elif script == 'effconn_ppi':

				config = PPI(subject=subject,
					     projectdir=config_args['projectdir'] 
					     processed_dir=config_args['processed_dir'] 
					     parent_dir=config_args['parent_dir']
					     stats_dir=config_args['stats_dir']
					     roi_type=config_args['roi_type']
					     roicenter_value=config_args['roicenter_value']
	        			     roiradius_value=config_args['roiradius_value'] 
					     session_value=config_args['session_value'] 
					     contrast_type=config_args['contrast_type'] 
					     contrastweights_value=config_args['contrastweights_value'])

			elif script == 'effconn_ppi_volrep'

				config = PPIVolRep(subject=subject,
						   projectdir=config_args['projectdir']
						   processed_dir=config_args['processed_dir']
	        				   pipeline_type=config_args['pipeline_type']
	        				   stats_dir=config_args['stats_dir']
	        				   roi_type=config_args['roi_type']
	        				   roicenter_value=config_args['roicenter_value'] 
	        				   roiradius_value=config_args['roiradius_value'] 
	        				   session_value=config_args['session_value'] 
	        				   contrast_type=config_args['contrast_type']
	        				   contrastweights_value=config_args['contrastweights_value'])

			elif script == 'effconn_phipi'

				config = PHIPI(subject=subject,
					       projectdir=config_args['projectdir'],
	        			       processed_dir=config_args['processed_dir'],
	        			       parent_dir=config_args['parent_dir'],
	        			       stats_dir=config_args['stats_dir'],
	        			       roiname_list=config_args['roiname_list'],
	        			       roicenter_list=config_args['roicenter_list'], 
	        			       roiradius_list=config_args['roiradius_list'],
	        			       session_value=config_args['session_value'],
	        			       contrast_type=config_args['contrast_type']
	        			       contrastweights_value=config_args['contrastweights_value'])

			elif script == 'scsnl_gPPI':

				config = gPPI(subject=subject,
					      project_dir=config_args['project_dir'],
					      parent_dir=config_args['parent_dir'],
					      stats_dir=config_args['stats_dir'],
					      pipeline=config_args['pipeline'],
					      roifile_list=config_args['roifile_list'],
					      roiname_list=config_args['roiname_list'],
					      tasks_to_include=config_args['tasks_to_include'],
					      confound_names=config_args['confound_names'])

			elif script == 'contrastchange':
	
				config = ContrastChange(subject=subject,
							projectdir=config_args['projectdir'],
							contrastfile=config_args['contrastfile'],
							statsdir=config_args['statsdir'],
							template_dir=config_args['template_dir'])

			input_args = config.__dict__
			savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
			submit_job(config_path,script,nameID,submit_script,slurm_args)








