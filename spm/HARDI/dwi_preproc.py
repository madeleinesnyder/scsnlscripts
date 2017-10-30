from nipype.interfaces import ants
from nipype.interfaces import fsl
from nipype.interfaces import freesurfer as fs
from nipype.interfaces import utility as niu
from nipype.interfaces import io as nio
from nipype.interfaces import mrtrix3 as mrt3
from nipype.interfaces.utility import Function
from nipype.pipeline import engine as pe
from nipype.interfaces import dipy
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
import os.path as op
import os

def dwidenoise(in_file, nthreads=1):
    from subprocess import call
    import os.path as op
    import os
    out_dir = os.getcwd()
    in_fname = op.basename(in_file)
    out_file = op.join(out_dir, in_fname.split('.')[0] + '_den.nii.gz') 
    print(out_file)
    call(['dwidenoise', '-nthreads', str(nthreads), in_file, out_file])
    return out_file

def dwiextract_b0(in_file, bval, bvec):
    from subprocess import call
    import os.path as op
    import os
    out_dir = os.getcwd()
    in_fname = op.basename(in_file)
    out_file = op.join(out_dir, in_fname.split('.')[0] + '_b0.nii.gz')
    #call(['dwiextract', '-bzero', '-grad', grad_file, in_file, out_file])
    call(['dwiextract', '-bzero', '-fslgrad', bvec, bval, in_file, out_file])
    return out_file

def dilate_nii(in_file,  npass=5):
    import numpy as np
    import nibabel as nib
    from scipy.ndimage.morphology import binary_dilation
    import os.path as op
    import os
    out_dir = os.getcwd()
    #out_dir = op.dirname(in_file)
    in_fname = op.basename(in_file)
    out_file = op.join(out_dir, in_fname.split('.')[0] + '_dil%i.nii.gz' % npass)
    #print(out_file)
    in_img = nib.load(in_file)
    data = in_img.get_data().astype(np.bool)
    data_dil = binary_dilation(data, iterations=npass)
    nib.save(nib.Nifti1Image(data_dil.astype(np.int16), in_img.get_affine()), out_file)
    return out_file

def prep_ecc_files(sidecar_file, dwi_file, in_bvec_file):
    #but c.f.: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;92136ade.1506
    # "Otherwise you can typically just go with 0.05. And also make sure that you use the same values for topup and eddy if you are using them together."
    import json, os 
    import os.path as op
    import nibabel as nib
    import numpy as np
    with open(sidecar_file, 'r') as f:
        sidecar = json.load(f)
    acqp_file = op.join(os.getcwd(),'acqp')
    acqp = sidecar['acqp']
    with open(acqp_file, 'w') as f:
        f.write(' '.join([str(x) for x in acqp]))
    index_file = op.join(os.getcwd(), 'index')
    nidx = nib.load(dwi_file).shape[-1]
    with open(index_file, 'w') as f:
        f.write(' '.join(nidx * ['1']))
    # Prep the bvecs:
    # see: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/FAQ#What_conventions_do_the_bvecs_use.3F
    # see: http://www.mrtrix.org/2016/04/15/bug-in-fsl-bvecs-handling/
    # FSL prefers LAS. Our bvecs are RAS. So we flip the x dim. 
    bvec_fname = os.path.basename(in_bvec_file)
    fsl_bvec_file = os.path.join(os.getcwd(), '%s_FSL' % bvec_fname)
    bvecs_fsl = np.loadtxt(in_bvec_file).T
    bvecs_fsl[0,:] *= -1
    np.savetxt(fsl_bvec_file, bvecs_fsl)
    return acqp_file, index_file, fsl_bvec_file

def dwibiascorr(in_file, grad_file):
    import os
    from subprocess import call
    in_fname = os.path.basename(in_file)
    out_dir = os.getcwd()
    out_file = os.path.join(out_dir, in_fname.split('.')[0] + '_dwibc.nii.gz')
    call(['dwibiascorrect', '-grad', grad_file, '-ants', in_file, out_file])
    return out_file

def bvecs2grads(in_bvec_file, in_bval_file, fsl_in = True):
    import numpy as np
    import os
    bvecs = np.loadtxt(in_bvec_file)
    bvals = np.loadtxt(in_bval_file)[:,np.newaxis]
    if fsl_in:
        bvecs = bvecs.T
        bvecs[:,0] *= -1
    grad_file = os.path.join(os.getcwd(), 'grads')
    np.savetxt(grad_file, np.c_[bvecs, bvals])
    return grad_file

def dwi_QA_calc(dwi_raw, eddy_params, grad_file, T=1000, D=0.001):
    import numpy as np
    import nibabel as nib
    import os
    motpars = np.loadtxt(eddy_params)[:,:6] #first 6 columns are motion correction params
    #1. Motion-based measures
    # Sum_k[sqrt(dx_k^2 + dy_k^2 + dz_k^2)] / Nvol
    avg_trans = np.r_[0,(np.diff(motpars[:,:3], axis=0) ** 2).sum(axis=1)].mean()
    # Sum_k[(|da_k| + |dg_k| + |db_k|)] / Nvol
    avg_rot = np.r_[0,np.sum(np.abs(motpars[:,3:6]),axis=1)].mean()
    #2. Signal-based measures 
    grads  = np.loadtxt(grad_file)
    bvals = grads[:,-1]
    bvecs = grads[:,:3]
    b = bvals.max()
    Imin = T * np.exp(-b * D)
    raw_img = nib.load(dwi_raw)
    dwidx = np.where(bvals > 1)[0]
    data = raw_img.get_data()
    # number of voxels with intensity > Imin in all DW volumes
    ri = np.squeeze(np.apply_over_axes(np.sum,(data[...,dwidx] > Imin).astype(np.float64), axes = [0,1]))
    # number of voxels with intensity > Imin in 1st DW volume
    r1 = np.squeeze(np.apply_over_axes(np.sum, (data[...,dwidx[0]] > Imin).astype(np.float64), axes=[0,1]))
    # Benner score
    Si = 2 - ri / (0.7 * r1[:,np.newaxis]) # slice x vol
    Si[ri < 0.05*np.prod(data.shape[:2])] = 0.0 #ignore slices where ri < 5% # voxels in the slice
    pct_bad = (Si > 1).astype(np.float64).sum() / np.prod(Si.shape) # percent bad slices
    if np.any(Si > 1):
        avg_drop = Si[Si>1].mean() # average dropout score
    else:
        avg_drop = 0.0
    qa_file = os.path.join(os.getcwd(), 'qa_measures')
    np.savetxt(qa_file, np.r_[avg_trans,avg_rot,pct_bad,avg_drop])
    return qa_file

## Maddie comments
## Main assemblor of all diffusion workflow steps 
#

def init_dwi_wf(opts, name='dwi_preproc_wf'):
    input_node = pe.Node(
        niu.IdentityInterface(fields=['t1', 'dwi', 'bval', 'bvec', 'bids_sidecar']),
        name='inputnode')
    input_node.inputs.dwi = opts.input_dwi
    input_node.inputs.bval = opts.bval
    input_node.inputs.bvec = opts.bvec
    input_node.inputs.bids_sidecar = opts.bids_sidecar
    input_node.inputs.t1 = opts.input_t1

    output_node = pe.Node(niu.IdentityInterface(
        fields=['dwi_mask','b0','grads','dt','fa','response','fODF', 'qa_file', 'dwi_processed',
                'eddy_params', 'xfm_anat2dwi', 'xfm_dwi2anat', 't1_processed']),
        name='outputnode')

    # 1. DENOISE
    dwi_den = pe.Node(Function(input_names=['in_file','nthreads'],
                               output_names=['out_file'],
                               function=dwidenoise),
                      name='dwi_denoise')
    dwi_den.inputs.nthreads = opts.nthreads

    # 2. GENERATE DWI BRAIN MASK
    # # Extract mean B0 (raw+denoised)
    dwi_get_b0 = pe.Node(Function(input_names=['in_file', 'bval', 'bvec'],
                                  output_names=['out_file'],
                                  function=dwiextract_b0),
                         name='dwiraw_ext_b0')

    # # Generate DWI space mask using FSL's bet2
    dwi_b0_mask = pe.Node(fsl.BET(frac = 0.35, mask=True, no_output=True),
                          name = 'dwi_b0_mask')

    # # Dilate the mask generously ( helpful for eddy [motion], no bad side effects at this stage )
    dwi_mask_dil = pe.Node(Function(input_names = ['in_file'],
                                    output_names=['out_file'],
                                    function=dilate_nii),
                           name='dwi_mask_dil')

    # Prep thea acqp and index files for Eddy
    prep_ecc = pe.Node(Function(input_names = ['sidecar_file', 'dwi_file', 'in_bvec_file'],
                                output_names = ['acqp_file', 'index_file', 'fsl_bvec_file'],
                                function=prep_ecc_files),
                       name = 'prep_data_ecc')

    # EDDY
    dwi_ecc = pe.Node(fsl.Eddy(num_threads = opts.nthreads, 
                               repol=True, 
                               environ = {'OMP_NUM_THREADS' : str(opts.nthreads),
                                          'FSLOUTPUTTYPE' : 'NIFTI_GZ'}), 
                      name = 'eddy')

    # Get flipped FSL-format bvecs and put them into mrtrix-style gradient table
    prep_grads = pe.Node(Function(input_names = ['in_bvec_file', 'in_bval_file'],
                                  output_names = ['grad_file'],
                                  function = bvecs2grads),
                         name = 'bvbv2grads')

    dwi_QA = pe.Node(Function(input_names = ['dwi_raw', 'eddy_params', 'grad_file'],
                              output_names = ['qa_file'],
                              function = dwi_QA_calc),
                     name = 'dwi_QA_calc')

    dwi_biascorr = pe.Node(Function(input_names  = ['in_file', 'grad_file'],
                                    output_names = ['out_file'],
                                    function = dwibiascorr),
                           name='dwi_bias_corr')

    dwi_reslice = pe.Node(dipy.Resample(vox_size = (2.,2.,2.), 
                                        interp=5), 
                          name='reslice_dwi')

    # Extract B0 from eddy+motion+bias corrected image
    dwi_corr_get_b0 = pe.Node(Function(input_names=['in_file', 'bval', 'bvec'],
                                      output_names=['out_file'],
                                      function=dwiextract_b0),
                             name='dwicorr_ext_b0')

    # Get the B0 brain image 
    dwi_corr_b0_mask = pe.Node(fsl.BET(frac=0.35, mask=True),
                              name = 'dwicorr_b0_mask')

    # Compute tensors
    dwi2dt = pe.Node(mrt3.FitTensor(nthreads = opts.nthreads, out_file = 'dt.nii.gz'),
                     name = 'dwi2tensor')

    dt2fa = pe.Node(mrt3.TensorMetrics(out_fa = 'fa.nii.gz'),
                    name = 'tensor2fa')

    # Setup the anat WF for (minimal) T1 processing and T1<->DWI Registration
    anat_wf = init_anat_wf(opts)

    # Estimate Response Function using 'Tournier' Algo
    dwi_est_RF = pe.Node(mrt3.ResponseSD(nthreads = opts.nthreads,
                                         max_sh = opts.lmax,
                                         args = 'tournier'
                                        ),
                         name = 'estimate_RF')
    
    # Estimate fODFs using CSD
    dwi_est_FOD = pe.Node(mrt3.EstimateFOD(nthreads = opts.nthreads,
                                           max_sh = opts.lmax,
                                           args = 'csd',
                                           out_file = 'fODF%i.nii.gz' % opts.lmax),
                          name = 'estimate_FOD_CSD')

    datasink = pe.Node(nio.DataSink(), name='sinker')
    datasink.inputs.base_directory = opts.output_basedir
    #datasink.inputs.container = op.join(opts.output_basedir, 'prep')

    workflow = pe.Workflow(name=name)
    workflow.connect([
        (input_node, dwi_den, [('dwi', 'in_file')]),
        (dwi_den, dwi_get_b0, [('out_file', 'in_file')]),
        (input_node, dwi_get_b0, [('bval','bval'), ('bvec','bvec')]),
        (dwi_get_b0, dwi_b0_mask, [('out_file', 'in_file')]),
        (dwi_b0_mask, dwi_mask_dil, [('mask_file', 'in_file')]),
        (dwi_den, output_node, [('out_file', 'dwi_den')]),
        (input_node, prep_ecc, [('bids_sidecar', 'sidecar_file'), 
                                ('dwi', 'dwi_file'),
                                ('bvec', 'in_bvec_file')]),
        (prep_ecc, dwi_ecc, [('acqp_file', 'in_acqp'), 
                             ('index_file', 'in_index'),
                             ('fsl_bvec_file', 'in_bvec')]),
        (dwi_den, dwi_ecc, [('out_file', 'in_file')]),
        (input_node, dwi_ecc, [('bval', 'in_bval')]),
        (dwi_mask_dil, dwi_ecc, [('out_file', 'in_mask')]),
        (dwi_ecc, dwi_biascorr, [('out_corrected', 'in_file')]),
        (dwi_ecc, prep_grads, [('out_rotated_bvecs', 'in_bvec_file')]), #TODO: test this
        (input_node, prep_grads, [('bval', 'in_bval_file')]),
        (dwi_ecc, dwi_QA, [('out_parameter', 'eddy_params')]),
        (dwi_ecc, output_node, [('out_parameter', 'eddy_params')]),
        (input_node, dwi_QA, [('dwi', 'dwi_raw')]),
        (prep_grads, dwi_QA, [('grad_file', 'grad_file')]),
        (dwi_QA, output_node, [('qa_file', 'qa_file')]),
        (prep_grads, output_node, [('grad_file', 'grads')]),
        (prep_grads, dwi_biascorr, [('grad_file', 'grad_file')]),
        (dwi_biascorr, dwi_reslice, [('out_file', 'in_file')]),        
        (dwi_reslice, dwi_corr_get_b0, [('out_file', 'in_file')]),
        (dwi_reslice, output_node, [('out_file', 'dwi_processed')]),
        (input_node, dwi_corr_get_b0, [('bval','bval'), ('bvec','bvec')]), # just needs the b0 indices
        (dwi_corr_get_b0, dwi_corr_b0_mask, [('out_file', 'in_file')]),
        (dwi_corr_b0_mask, dwi2dt, [('mask_file', 'in_mask')]),
        (dwi_corr_b0_mask, output_node, [('out_file', 'b0')]),
        (dwi_reslice, dwi2dt, [('out_file', 'in_file')]),
        (prep_grads, dwi2dt, [('grad_file', 'grad_file')]),
        (dwi2dt, output_node, [('out_file', 'dt')]),
        (dwi2dt, dt2fa, [('out_file','in_file')]),
        (dwi_corr_b0_mask, dt2fa, [('mask_file', 'in_mask')]),
        (dt2fa, output_node, [('out_fa', 'fa')]),
        (dwi_reslice, dwi_est_RF, [('out_file', 'in_file')]), 
        (prep_grads, dwi_est_RF, [('grad_file', 'grad_file')]),
        (dwi_corr_b0_mask, dwi_est_RF, [('mask_file', 'in_mask')]),
        (dwi_est_RF, output_node, [('out_file',  'SD_RF')]),
        (dwi_est_RF, dwi_est_FOD, [('out_file', 'response')]),
        (dwi_est_RF, output_node, [('out_file', 'response')]),
        (dwi_reslice, dwi_est_FOD, [('out_file', 'in_file')]),
        (dwi_corr_b0_mask, dwi_est_FOD, [('mask_file', 'in_mask')]),
        (prep_grads, dwi_est_FOD, [('grad_file', 'grad_file')]),
        (dwi_est_FOD, output_node, [('out_file', 'fODF')]),
        (dt2fa, anat_wf, [('out_fa', 'inputnode.fa')]),
        (input_node, anat_wf, [('t1','inputnode.t1')]),
        (anat_wf, output_node, [('outputnode.xfm_anat2dwi', 'xfm_anat2dwi'),
                                ('outputnode.xfm_dwi2anat', 'xfm_dwi2anat'),
                                ('outputnode.t1_processed', 't1_processed')]),
        (output_node, datasink, [('dwi_mask', 'prep.@mask'),
                                 ('b0', 'prep.@b0'),
                                 ('dt', 'prep.@dt'),
                                 ('fa', 'prep.@fa'),
                                 ('fODF', 'prep.@FOD'),
                                 ('response', 'prep.@response'),
                                 ('eddy_params', 'prep.@eddy_params'),
                                 ('dwi_processed', 'prep.@dwi_processed'), 
                                 ('grads', 'prep.@grads'),
                                 ('qa_file', 'prep.@QA'),
                                 ('xfm_anat2dwi', 'prep.xfm.@anat2dwi'),
                                 ('xfm_dwi2anat', 'prep.xfm.@dwi2anat'),
                                 ('t1_processed', 'prep.@t1')])
        ])
    return workflow

def init_anat_wf(opts, name='anat_wf'):

    workflow = pe.Workflow(name=name)

    input_node = pe.Node(
        niu.IdentityInterface(fields=['t1', 'fa']),
        name='inputnode')

    output_node = pe.Node(niu.IdentityInterface(
        fields=['xfm_anat2dwi', 'xfm_dwi2anat', 't1_processed']),
        name='outputnode')

    skullstrip_wf = init_skullstrip_watershed_wf(n4_nthreads = opts.nthreads)

    t1_dwi_reg = pe.Node(ants.Registration(),
                         name = 't1_2_dwi')

    t1_dwi_reg.inputs.transforms = ['Rigid', 'Affine','SyN']
    t1_dwi_reg.inputs.interpolation = 'BSpline'
    t1_dwi_reg.inputs.transform_parameters = [(0.1,), (0.1,), (0.1, 3.0, 0.05)]
    t1_dwi_reg.inputs.number_of_iterations = [[1000, 500, 250, 100], [500, 250, 100], [100, 50, 30]]
    t1_dwi_reg.inputs.dimension = 3
    t1_dwi_reg.inputs.write_composite_transform = True
    t1_dwi_reg.inputs.collapse_output_transforms = True
    t1_dwi_reg.inputs.metric = ['MI']*3
    t1_dwi_reg.inputs.metric_weight = [1]*3 # only matters if using more than 1 metric per step
    t1_dwi_reg.inputs.radius_or_number_of_bins = [32]*3
    t1_dwi_reg.inputs.sampling_strategy = ['Regular']*3
    t1_dwi_reg.inputs.sampling_percentage = [0.25]*3
    t1_dwi_reg.inputs.convergence_threshold = [1e-6, 1.e-6, 1.e-8]
    t1_dwi_reg.inputs.convergence_window_size = [20, 20, 10]
    t1_dwi_reg.inputs.smoothing_sigmas = [[3,2,1,0], [2,1,0], [2,1,0]]
    t1_dwi_reg.inputs.sigma_units = ['vox'] * 3
    t1_dwi_reg.inputs.shrink_factors = [[8,4,2,1], [4,2,1], [3,2,1]]
    t1_dwi_reg.inputs.use_estimate_learning_rate_once = [True] * 3
    t1_dwi_reg.inputs.use_histogram_matching = [True] * 3
    t1_dwi_reg.inputs.num_threads = opts.nthreads
    t1_dwi_reg.inputs.output_warped_image = 'FA_T1-space.nii.gz'
    t1_dwi_reg.inputs.output_transform_prefix = "xfm_DW2T1_"

    workflow.connect([
        (input_node, skullstrip_wf, [('t1', 'inputnode.in_file')]),
        (skullstrip_wf, t1_dwi_reg, [('outputnode.out_file', 'fixed_image')]),
        (skullstrip_wf, output_node, [('outputnode.out_file', 't1_processed')]),
        (input_node, t1_dwi_reg, [('fa', 'moving_image')]),
        (t1_dwi_reg, output_node, [('composite_transform', 'xfm_dwi2anat'),
                                  ('inverse_composite_transform', 'xfm_anat2dwi')])
        ])
    return workflow


def init_skullstrip_watershed_wf(name='HWASkullStripWorkflow', n4_nthreads=1):
    """
    Skull-stripping workflow utilizing FreeSurfer's mri_watershed
    """

    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file']),
                        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['bias_corrected', 'out_file', 'out_mask', 'bias_image']), 
        name='outputnode')

    t1_denoise = pe.Node(ants.DenoiseImage(dimension = 3), name = 'denoise')

    n4_correct = pe.Node(
        ants.N4BiasFieldCorrection(dimension=3, save_bias=True, num_threads=n4_nthreads),
        name='n4')

    t1_skull_strip = pe.Node(fs.WatershedSkullStrip(out_file = 'T1w_corrected_brain.nii.gz'),
                             name='watershed')
    #t1_skull_strip.inputs.out_file = 'T1w_corrected_brain.nii.gz'

    create_mask = pe.Node(fs.Binarize(
        min=1.e-3, dilate=0, out_type='nii.gz', binary_file = 'T1w_corrected_brain_mask.nii.gz'),
        name='binarize')
    #create_mask.inputs.binary_file = 'T1w_corrected_brain_mask.nii.gz'

    workflow.connect([
        (inputnode, t1_denoise, [('in_file', 'input_image')]),
        (t1_denoise, n4_correct, [('output_image', 'input_image')]),
        (n4_correct, t1_skull_strip, [('output_image', 'in_file')]),
        (n4_correct, outputnode, [('output_image', 'bias_corrected'),
                                  ('bias_image', 'bias_image')]),
        (t1_skull_strip, create_mask, [('out_file', 'in_file')]),
        (create_mask, outputnode, [('binary_file', 'out_mask')]),
        (t1_skull_strip, outputnode, [('out_file', 'out_file')])
    ])
    return workflow

def get_parser():
    parser = ArgumentParser(description='SCSNL HARDI pipeline',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument('project_dir', action='store',
                        help='Path to SCSNL-style project directory')
    parser.add_argument('PID', action='store', type=int)
    parser.add_argument('visit', action='store', type=int)
    parser.add_argument('session', action='store', type=int)
    parser.add_argument('--lmax', dest = 'lmax', default = 8,
                        help = 'Maximum harmonic degree for CSD')
    parser.add_argument('--nthreads', dest = 'nthreads', default = 16,
                        help = 'Number of cores to use for multi-threaded stages')
    return parser

## TODO 
#argparse.... 
#opts.lmax
#opts.subject
#opts.bids_sidecar = get sidecar from subject dir
#opts.bval = get bval from subject dir
#opts.bvec = get bvec from subject dir
#opts.nthreads = int(os.environ['SLURM_NTASKS'])
#opts.output_basedir = #TODO output_basedir

def main():
    opts = get_parser().parse_args()
    opts.subject_dir = op.join(opts.project_dir, 'data', 'imaging', 'participants', str(opts.PID), 
        'visit%i' % opts.visit, 'session%i' % opts.session)
    opts.bids_sidecar = op.join(opts.subject_dir, 'dwi', 'dwi_raw.json')
    opts.bval = op.join(opts.subject_dir, 'dwi', 'dwi_raw.bval')
    opts.bvec = op.join(opts.subject_dir, 'dwi', 'dwi_raw.bvec')
    opts.input_dwi = op.join(opts.subject_dir, 'dwi', 'dwi_raw.nii.gz')
    opts.input_t1 = op.join(opts.subject_dir, 'anatomical', 'T1w.nii.gz')
    if not opts.nthreads:
        opts.nthreads = 16 #int(os.environ['SLURM_NTASKS'])
    opts.output_basedir = op.join(opts.subject_dir, 'dwi')
    workflow = init_dwi_wf(opts)

    workflow.write_graph(format='svg', dotfilename='diffusion_wf.dot')
    #workflow.run(plugin='SLURMGraph', plugin_args = {'dont_resubmit_completed_jobs': True, 'sbatch_args' : ' '})
    workflow.base_dir = op.join(opts.subject_dir, 'dwi', 'wf')
    #'/home/users/kochalka/hardi/diffusion_wf'
    workflow.run()

if __name__ == '__main__':
    main()

# TODO: 
#if reg_method == 'bbr' and opts.subject_freesurfer_dir:
#
#elif:  
#

"""
t1_dwi = pe.Node(
    ANTSRegistrationRPT(
        generate_report=True,
        num_threads=omp_nthreads,
        from_file='/path/to/t1-fa_registration_precise_000.json' #TODO
        #flavor='testing' if debug else 'precise',
    ),
    name='t1_2_dwi'
)
"""

# TODO: finalize reslicing... 


