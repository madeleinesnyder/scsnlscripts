import os
import itertools
import nibabel as nib
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def dice(im1, im2, empty_score=1.0):
    """
    Computes the Dice coefficient, a measure of set similarity.
    Parameters
    ----------
    im1 : array-like, bool
        Any array of arbitrary size. If not boolean, will be converted.
    im2 : array-like, bool
        Any other array of identical size. If not boolean, will be converted.
    Returns
    -------
    dice : float
        Dice coefficient as a float on range [0,1].
        Maximum similarity = 1
        No similarity = 0
        Both are empty (sum eq to zero) = empty_score
        
    Notes
    -----
    The order of inputs for `dice` is irrelevant. The result will be
    identical if `im1` and `im2` are switched.
    """
    im1 = np.asarray(im1).astype(np.bool)
    im2 = np.asarray(im2).astype(np.bool)

    if im1.shape != im2.shape:
        raise ValueError("Shape mismatch: im1 and im2 must have the same shape.")

    im_sum = im1.sum() + im2.sum()
    if im_sum == 0:
        return empty_score

    # Compute Dice coefficient
    intersection = np.logical_and(im1, im2)

    return 2. * intersection.sum() / im_sum

img_dir = '/oak/stanford/groups/menon/scsnlscripts/ica-aroma/analysis/mathFUN_groupstats'
timepoints = ['pre','post']
grps = ['ALL_no_mvmnt_criteria','10mm_max_displacement','6mm_max_displacement','3mm_max_displacement']
contrasts = ['001T_1_Dot_near-far','003T_3_Dot_all_correct-rest','005T_5_Num_near-far','007T_7_Num_all_correct-rest']


f, axes = plt.subplots(2,4)
f.suptitle('Dice index for mathFUN group stats',fontsize=16)
for time_i, time in enumerate(timepoints):
    for con_i, con in enumerate(contrasts):
        imgdata = []
        for grp in grps:

            negfile = time+'_'+grp+'_'+'p01_128voxels_'+con+'_neg_p_unc.01_128vox.nii'
            posfile = time+'_'+grp+'_'+'p01_128voxels_'+con+'_pos_p_unc.01_128vox.nii'
            negfile = os.path.join(img_dir,time,grp,con,negfile)
            posfile = os.path.join(img_dir,time,grp,con,posfile)

            posimg = nib.load(posfile)
            posdat = posimg.get_data()
            imgdata.append(posdat)

            dicemat = np.zeros([len(imgdata),len(imgdata)])
            for (i1,l1), (i2,l2) in itertools.combinations(enumerate(imgdata), 2):
                dicemat[i1,i2] = dice(l1,l2)
                dicemat[i2,i1] = dice(l1,l2)

        if con_i == 3:
            showcbar = True
            cbar_ax = f.add_axes([0.91,0.3,0.03,0.4])
        else:
            showcbar = False
            cbar_ax = None
        hm = sns.heatmap(dicemat,xticklabels=['All','10mm','6mm','3mm'],yticklabels=['All','10mm','6mm','3mm'],
                         ax=axes[time_i,con_i],vmin=0,vmax=1,cbar=showcbar,cbar_ax=cbar_ax,square=True)
        axes[time_i,con_i].set_title(con)
axes[0,0].set_ylabel('Pre')
axes[1,0].set_ylabel('Post')

plt.show()
















