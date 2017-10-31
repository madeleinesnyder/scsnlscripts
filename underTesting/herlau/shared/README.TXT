VERSION 0.1

MDS-MCMC Instructions:

INSTALLING:
Simply unzip the main zip file. The code is in pure matlab except two external dependencies which should be compiled.
First build the mmx 3d matrix multiplication library by running 

code/GPMDS/mmx/build_mmx.m

In case you wish to run the supplied spm-DCM wrapper you also need to build the included version of SPM12 located in /spm12. 
Use the standard online instructions to build SPM12. 

TESTING:
see code/GPMDS_example/GMPDS_Example*.m for examples of running the method.
The examples show how the MDS-MCMC method can be evaluated on fMRI data (in this case, the opto-rat dataset),
but also run DCM and MDS. Output format should be reasonably self-explanatory and I believe the wrapper is correct; however 
use at own peril. 

SPECIFYING YOUR OWN EXPERIMENTS:
The example script should give an idea for the standard data format. The main line is:

[stats,p2] = GPMDS_matlab(p,opts)

where p is the data structure for the fMRI data and opts are options for the sampler (mainly number of iterations). 
As output, stats are (summary) statistics recorded by the sampler and p2 is the (current) state of the chain (to support resuming the sampling process). 

To provide a brief overview about the data structure of p:

p.C 	<- M x M x J interaction matrix for M regions over J conditions
p.s(n) 	<- data structure for subject n.
p.dT 	<- temporal resolution (TR) of the fMRI experiment

Let's examine this in more details:

p.s(n).y <- m x T observed fMRI data
p.s(n).v <- J x T binary matrix indicating which of the J C-matrices are active. 
p.s(n).u <- M x T x S Entry u(m,:,s) is a binary vector indicating the presences of stimuli condition s. These can be independent of v. 
It is recommended to set p.s(n).u(:,:,1) = 1 to indicate there may be a constant base-rate activity different than zero and then specify additional conditions in the following blocks. For instance p.s(n).u(4, 50:100, 2) = 1 indicates that this subject was stimulated in region 4
during observation 50-100. 

It is recommended to run the sampler for a long time. A really, really long time. For that reason the main scripts are designed to allow interruption and resuming. Ask for scripts that should be cluster-ready. It is recommended to run the method several times to test for stability. 
