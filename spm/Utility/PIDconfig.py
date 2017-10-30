from subprocess import call
from scipy.io import savemat
import numpy as np
from numpy import array
import os
import sys

oakpath = os.environ['OAK']

class Config(object):
	'''Basic config class

	Attributes:
		subject: A subject to run your script on
		runlist: A list of runs to run your script on
		datadir: Directory holding raw data
	'''
	datadir = oakpath + 'rawdata'
	# self assign oakpath for bash template

	def __init__(self, subject,runlist):
		self.subject = subject
		self.runlist = runlist


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
		- runlist
		- pipeline
		- project directory
	  Optional:
		- inputimgprefix (defaults to None)
		- SPGRsubjectlist (defaults to None)
		- SPGRfilename (defaults to 'watershed_spgr')
		- tr_val (defaults to 2.0)
	'''

	def __init__(self,subject,runlist,pipeline,projectdir,
				 inputimgprefix='',SPGRsubjectlist='',
				 SPGRfilename='watershed_spgr',tr_val=2.0):
		super(Preprocessing, self).__init__(subject,runlist)
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

	def __init__(self,subject,runlist,pipeline_type,project_dir,taskdesign_file,
		     contrast_file,movement_include,fmri_type,stats_dir='stats_spm8',
		     preprocessed_dir='smoothed_spm8',parent_dir='',image_type='nii',
		     volrepair_include=0,volpipeline_type='',volrepaired_dir='',repairedstats_dir='',
		     template_dir = oakpath + '/scsnlscripts/spm/BatchTemplates'):
		super(IndividualStats, self).__init__(subject,runlist)
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
		if self.pipeline_type not in pipeline_options:
			raise Exception('Incorrect Pipeline specified')

		if not os.path.isdir(self.project_dir):
			raise Exception('Project directory does not exist')

class MultiIndividualStats(Config):
	'''MultiIndividualStats config class

		TODO: Logic for asserting subject lists are the same length
	'''

	def __init__(self,subject,runlist,pipeline_type,project_dir,taskdesign_file,
		     contrast_file,movement_include,fmri_type,stats_dir='stats_spm8',
		     preprocessed_dir='smoothed_spm8',image_type='nii',
		     volrepair_include=0,volpipeline_type='',volrepaired_dir='',repairedstats_dir='',
		     template_dir = oakpath +'scsnlscripts/spm/BatchTemplates'):
		super(MultiIndividualStats, self).__init__(subject,runlist)
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

	def __init__(self,subject,runlist,project_dir,pipeline,TR,
		     roi_dir,roi_list,fmri_type,stats_dir,preprocessed_dir,
		     numtrunc=[0,0],data_type='nii',nonyear_dir='',bandpass=0,fl=0,fh=0.2):
		super(FunctionalConnectivity, self).__init__(subject,runlist)
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
		     run_value, contrast_type, contrastweights_value):

		self.subject = subject
		self.projectdir = projectdir 
		self.processed_dir = processed_dir 
		self.parent_dir = parent_dir
		self.stats_dir = stats_dir
		self.roi_type = roi_type
		self.roicenter_value = roicenter_value
	        self.roiradius_value = roiradius_value 
		self.run_value = run_value 
		self.contrast_type = contrast_type 
		self.contrastweights_value = contrastweights_value

class PPIVolRep(object):
	'''PPIVolRep config class
	'''

	def __init__(self, projectdir, processed_dir, pipeline_type,
		     stats_dir, roi_type, roicenter_value, roiradius_value, 
		     run_value, contrast_type, contrastweights_value):

		self.subject = subject
		self.projectdir = projectdir
		self.processed_dir = processed_dir
	        self.pipeline_type = pipeline_type
	        self.stats_dir = stats_dir
	        self.roi_type = roi_type
	        self.roicenter_value = roicenter_value 
	        self.roiradius_value = roiradius_value 
	        self.run_value = run_value 
	        self.contrast_type = contrast_type
	        self.contrastweights_value = contrastweights_value

class PHIPI(object):
	'''PHIPI config class
	'''

	def __init__(self, subject, projectdir, processed_dir, parent_dir,
	             stats_dir, roiname_list, roicenter_list, roiradius_list,
	             run_value, contrast_type, contrastweights_value):
		
		self.subject = subject
		self.projectdir = projectdir
	        self.processed_dir = processed_dir
	        self.parent_dir = parent_dir
	        self.stats_dir = stats_dir
	        self.roiname_list = roiname_list
	        self.roicenter_list = roicenter_list 
	        self.roiradius_list = roiradius_list
	        self.run_value = run_value
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
		     template_dir= oakpath +'scsnlscripts/spm/BatchTemplates'):
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
	   - runlist
	   - edatfolder
	   - onsetcolumn
	   - offsetcolumn
	   - taskname
	   - matchtaskname
	   - taskdesign
	   - restexist
	'''

	def __init__(self,subject,runlist,edatfolder,
		     onsetcolumn,restexist,taskname,rawdir, projectdir,
		     matchtaskname,taskdesignfilename,offsetcolumn,
		     globalonsetcolumn,globalonsetshift,
		     durationcolumn,presetduration,tasktypecolumn,runname):
		self.subject = subject
		self.runlist = runlist
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
		self.runname = runname

class RSAIndividual(object):

	def __init__(self,subject,roi_folder,roi_list,stats_folder,
	             map_type,map_pair_index,projectdir,output_dir,user_fname,marsbar_path):
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
	        self.projectdir = projectdir
	        self.output_dir = output_dir
	        self.user_fname = user_fname
	        self.marsbar_path = marsbar_path


class RSAWholeBrain(object):
	
	def __init__(self,projectdir,subject,statsfolder,maptype,mapindex,maskfile,outputdir,searchshape,searchradius):
		self.projectdir = projectdir
	        self.subject = subject
	        self.statsfolder = statsfolder
	        self.maptype = maptype
	        self.mapindex = array(mapindex)
	        self.maskfile = maskfile
	        self.outputdir = outputdir
		self.searchshape = searchshape
		self.searchradius = searchradius

class RSAGroup(object):

	def __init__(self,group_path,roi_folder,roi_list,map_folder,map_type,
	             map_index,projectdir,output_dir,user_fname,marsbar_path):

		self.group_path = group_path
	        self.roi_folder = roi_folder
		if type(roi_list) == list:
	        	self.roi_list = array(roi_list)
		else:
			self.roi_list = roi_list
	        self.map_folder = map_folder
	        self.map_type = map_type
	        self.map_index = map_index
	        self.projectdir = projectdir
	        self.output_dir = output_dir
	        self.user_fname = user_fname
	        self.marsbar_path = marsbar_path


class ContrastChange(object):
	'''ContrastChange config class
	'''

	def __init__(self, subject, projectdir,
		     contrastfile, statsdir,
		     template_dir= oakpath + 'scsnlscripts/spm/BatchTemplates'):
		self.subject = subject
		self.projectdir = projectdir
		self.contrastfile = contrastfile
		self.statsdir = statsdir
		self.template_dir = template_dir

class ContrastGenerator(object):
	'''ContrastGenerator config class
	'''

	def __init__(self, projectdir, numcontrasts, numruns, comparewithin, 
                   movementcorrection, contrastnames, contrast):
		self.projectdir = projectdir
		self.numcontrasts = numcontrasts
		self.numruns = numruns
		self.comparewithin = comparewithin
		self.movementcorrection = movementcorrection
		self.contrastnames = array(contrastnames, dtype=np.object)
		self.contrast = array(contrast, dtype=np.object)

class MovementExclusion(object):
	'''MovementExclusion config class
	'''
	
	def __init__(self, subject, runlist, rawdir, projectdir,preprocessedfolder,
			scantoscancrit):
		self.subject = subject
		self.runlist = runlist
		self.rawdir = rawdir
		self.projectdir = projectdir
		self.preprocessedfolder = preprocessedfolder
		self.scantoscancrit = scantoscancrit

class ROISignalLevel(object):

	
	def __init__(self, subject, projectdir,parent_folder,roi_folder,roi_list,stats_folder,
		     tscore_threshold,roi_result_folder,marsbar_path):
		self.subject = subject,
	        self.parent_folder = array(parent_folder),
	        self.roi_folder = roi_folder,
	        self.roi_list = roi_list,
	        self.stats_folder = stats_folder,
	        self.tscore_threshold = tscore_threshold,
	        self.projectdir = projectdir,
	        self.roi_result_folder = roi_result_folder,
	        self.marsbar_path = marsbar_path

class DicomToSpgr(object):
	'''DicomToSpgr config class
	'''

	def __init__(self, subject, projectdir, rawdir, numdicom):
		self.subject = subject
		self.projectdir = projectdir
		self.rawdir = rawdir
		self.numdicom = numdicom

class DTIACPC(object):

	def __init__(self,subject,projectdir,rawdir,parent_folder,participant_path,type_folder,vista_path):
		self.subject = subject
		self.projectdir = projectdir
		self.rawdir = rawdir
		self.parent_folder = parent_folder
		self.participant_path = participant_path
		self.type_folder = type_folder
		self.vista_path = vista_path

class DTIDCM2NII(object):

	def __init__(self,subject,projectdir,rawdir,parent_folder,numdicom,participant_path,type_folder,vist_path):
		self.subject = subject
		self.projectdir = projectdir
		self.rawdir = rawdir
		self.parent_folder = parent_folder
		self.numdicom = numdicom
		self.participant_path = participant_path
		self.type_folder = type_folder
		self.vista_path = vista_path

class DTIPreprocess(object):

	def __init__(self,subject,projectdir,participant_path,outname,vista_path,outBaseName,rawprepreocess_flag,
	             rawfittensor_flag,clobber,assetflag,numBootStrapSamples,eddyCorrect,excludeVols,bsplineInterpFlag):
		self.subject = subject,
		self.projectdir = projectdir
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
	
	def __init__(self,projectdir,rawdir,parent_folder,subject,vista_path,
	             type_folder,roilist,output_folder):
		self.projectdir = projectdir
	        self.rawdir = rawdir
	        self.parent_folder = parent_folder
	        self.subject = subject
		self.vista_path = vista_path
	        self.roilist = roi_list
	        self.type_folder = type_folder
	        self.output_folder = output_folder

class DTINiftiROI2DTIROI(object):

	def __init__(self,projectdir,rawdir,parent_folder,subject,vista_path,rois,type_folder):
		self.projectdir = projectdir
	        self.rawdir = rawdir
	        self.parent_folder = parent_folder
	        self.subject = subject
	        self.vista_path = vista_path
	        self.rois = rois
	        self.type_folder = type_folder

class DTITrackStats(object):

	def __init__(self,subject,projectdir,rawdir,parent_folder,subjectlist,
		     vista_path,roilist,type_folder,output_folder):
		self.subject = subject
		self.projectdir = projectdir
	        self.rawdir = rawdir
	        self.parent_folder = parent_folder
	        self.vista_path = vista_path
	        self.roilist = roilist
	        self.type_folder = type_folder
	        self.output_folder = output_folder

#class dti_wholebraintrack(object):

#	def __init__(self,subject,projectdir,rawdir,parent_folder,vista_path,type_folder,
#	             preprocessed_folder,faThresh,stepSizeMm,faThreshopts,lengthThreshMm,
#	             angleThresh,wPuncture,whichAlgorithm,whichInterp,seedVoxelOffsets):


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
	submit_script = '/oak/stanford/groups/menon/scsnlscripts/spm/Utility/submit_job.sh'

	call('mkdir %s'%(config_path),shell=True)

	if script == 'rsa_group':
		config = RSAGroup(group_path = config_args['group_path'],
				  roi_folder = config_args['roi_folder'],
				  roi_list = config_args['roi_list'],
				  map_folder = config_args['map_folder'],
				  map_type = config_args['map_type'],
	             		  map_index = config_args['map_index'],
				  projectdir = config_args['project_dir'],
				  output_dir = config_args['output_dir'],
				  user_fname = config_args['user_fname'],
				  marsbar_path = config_args['marsbar_path'])
		nameID = ''
		input_args=config.__dict__
		savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
		submit_job(config_path,script,nameID,submit_script,slurm_args)


	if script == 'contrast_generator':
		config = ContrastGenerator(projectdir = config_args['projectdir'],
					   numcontrasts = config_args['numcontrasts'],
					   numruns = config_args['numruns'],
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
						       runlist=config_args['runlist'],
						       tr_val=config_args['tr_val'],
						       inputimgprefix=config_args['inputimgprefix'],
						       pipeline=config_args['pipeline'],
						       SPGRfilename=config_args['SPGRfilename'],
						       projectdir=config_args['projectdir'])
				config.set_datadir(config_args['rawdir'])

			elif script == 'individualstats':

				config = IndividualStats(subject=subject,
							 image_type=config_args['image_type'],
							 pipeline_type=config_args['pipeline_type'],
							 project_dir=config_args['project_dir'],
							 parent_dir=config_args['parent_dir'],
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
				config.set_datadir(config_args['rawdir'])

			elif script == 'multi_individualstats':

				config = MultiIndividualStats(subject=subject,
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
				config.set_datadir(config_args['rawdir'])


			elif script == 'movement_exclusion':

				config = MovementExclusion(subject=subject,
							    runlist=config_args['runlist'],
							    rawdir=config_args['rawdir'],
							    projectdir=config_args['projectdir'],
							    preprocessedfolder=config_args['preprocessedfolder'],
							    scantoscancrit=config_args['scantoscancrit'])

			elif script == 'taskdesign_generation':

				config = TaskDesign(subject=subject,
						    runlist=config_args['runlist'],
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
						    runname=config_args['runname'],
						    taskdesignfilename=config_args['taskdesignfilename'],
						    restexist=config_args['restexist'])
			
			elif script == 'fconnect':

				config = FunctionalConnectivity(subject=subject,
							        data_type=config_args['data_type'],
							        nonyear_dir=config_args['nonyear_dir'],
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
							        roi_list=config_args['roi_list'],
							        numtrunc=config_args['numtrunc'],
							        fmri_type=config_args['fmri_type'])
				config.set_datadir(config_args['rawdir'])

			elif script == 'effconn_ppi':

				config = PPI(subject=subject,
					     projectdir=config_args['projectdir'], 
					     processed_dir=config_args['processed_dir'], 
					     parent_dir=config_args['parent_dir'],
					     stats_dir=config_args['stats_dir'],
					     roi_type=config_args['roi_type'],
					     roicenter_value=config_args['roicenter_value'],
	        			     roiradius_value=config_args['roiradius_value'],
					     run_value=config_args['run_value'], 
					     contrast_type=config_args['contrast_type'], 
					     contrastweights_value=config_args['contrastweights_value'])

			elif script == 'effconn_ppi_volrep':

				config = PPIVolRep(subject=subject,
						   projectdir=config_args['projectdir'],
						   processed_dir=config_args['processed_dir'],
	        				   pipeline_type=config_args['pipeline_type'],
	        				   stats_dir=config_args['stats_dir'],
	        				   roi_type=config_args['roi_type'],
	        				   roicenter_value=config_args['roicenter_value'], 
	        				   roiradius_value=config_args['roiradius_value'],
	        				   run_value=config_args['run_value'], 
	        				   contrast_type=config_args['contrast_type'],
	        				   contrastweights_value=config_args['contrastweights_value'])

			elif script == 'effconn_phipi':

				config = PHIPI(subject=subject,
					       projectdir=config_args['projectdir'],
	        			       processed_dir=config_args['processed_dir'],
	        			       parent_dir=config_args['parent_dir'],
	        			       stats_dir=config_args['stats_dir'],
	        			       roiname_list=config_args['roiname_list'],
	        			       roicenter_list=config_args['roicenter_list'], 
	        			       roiradius_list=config_args['roiradius_list'],
	        			       run_value=config_args['run_value'],
	        			       contrast_type=config_args['contrast_type'],
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

			elif script == 'rsa_individual':

				config = RSAIndividual(subject=subject,
						       roi_folder=config_args['roi_folder'],
						       roi_list=config_args['roi_list'],
						       stats_folder=config_args['stats_folder'],
	             				       map_type=config_args['map_type'],
						       map_pair_index=config_args['map_pair_idex'],
						       projectdir=config_args['projectsdir'],
						       output_dir=config_args['output_dir'],
						       user_fname=config_args['user_fname'],
						       marsbar_path=config_args['marsbar_path'])


			elif script == 'rsa_wholebrain':

				config = RSAWholeBrain(subject=subject,
						       projectdir=config_args['projectdir'],
						       statsfolder=config_args['statsfolder'],
						       maptype=config_args['maptype'],
						       mapindex=config_args['mapindex'],
						       maskfile=config_args['maskfile'],
						       outputdir=config_args['outputdir'],
						       searchshape=config_args['searchshape'],
						       searchradius=config_args['searchradius'])

			elif script == 'roi_signallevel':

				config = ROISignalLevel(subject=subject,
	        					parent_folder=config_args['parent_folder'],
	        					roi_folder=config_args['roi_folder'],
	        					roi_list=config_args['roi_list'],
	        					stats_folder=config_args['stats_folder'],
	        					tscore_threshold=config_args['tscore_threshold'],
	        					projectdir=config_args['projectdir'],
	        					roi_result_folder=config_args['roi_result_folder'],
	        					marsbar_path=config_args['marsbar_path'])

			elif script == 'contrastchange':
	
				config = ContrastChange(subject=subject,
							projectdir=config_args['projectdir'],
							contrastfile=config_args['contrastfile'],
							statsdir=config_args['statsdir'],
							template_dir=config_args['template_dir'])

			elif script == 'dicomtospgr':

				config = DicomToSpgr(subject=subject,
						     projectdir=config_args['projectdir'],
						     rawdir=config_args['rawdir'],
						     numdicom=config_args['numdicom'])

			elif script == 'dti_acpc':

				config = DTIACPC(subject=subject,
						 projectdir=config_args['projectdir'],
						 rawdir=config_args['rawdir'],
						 parent_folder=config_args['parent_folder'],
						 participant_path=config_args['participant_path'],
						 type_folder=config_args['type_folder'],
						 vista_path=config_args['vista_path'])

			elif script == 'dti_dcm2nii':

				config = DTIDCM2NII(subject=subject,
						    projectdir=config_args['projectdir'],
						    rawdir=config_args['rawdir'],
						    parent_folder=config_args['parent_folder'],
						    numdicom=config_args['numdicom'],
						    participant_path=config_args['participant_path'],
						    type_folder=config_args['type_folder'],
						    vista_path=config_args['vista_path'])

			elif script == 'dti_preprocess':

				config = DTIPreprocess(subject=subject,
						       projectdir=config_args['projectdir'],
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

				config = DTIBetweenROITrack(subject=subject,
						            projectdir=config_args['projectdir'],
	       						    rawdir=config_args['rawdir'],
	       						    parent_folder=config_args['parent_folder'],
	       						    vista_path=config_args['vista_path'],
	       						    roilist=config_args['type_folder'],
	       						    type_folder=config_args['type_folder'],
	       						    output_folder=config_args['output_folder'])

			elif script == 'dti_niftiroi2dtiroi':

				config = DTINiftiROI2DTIROI(subject=subject,
					                    projectdir=config_args['projectdir'],
	       						    rawdir=config_args['rawdir'],
	       						    parent_folder=config_args['parent_folder'],
	       						    vista_path=config_args['vista_path'],
	       						    rois=config_args['rois'],
	       						    type_folder=config_args['type_folder'])

			elif script == 'dti_trackstat':

				config = DTITrackStats(subject=subject,
						       projectdir=config_args['projectdir'],
	       					       rawdir=config_args['rawdir'],
	       					       parent_folder=config_args['parent_folder'],
	       					       vista_path=config_args['vista_path'],
	       					       roilist=config_args['roilist'],
	       					       type_folder=config_args['type_folder'],
	       					       output_folder=config_args['output_folder'])
						    

			input_args = config.__dict__
			savemat('%s/config_args_%s_%s.mat'%(config_path,script,nameID),input_args)
			submit_job(config_path,script,nameID,submit_script,slurm_args)








