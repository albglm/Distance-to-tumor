#!/bin/bash

#################################################################################################################
## Single patient / across timepoints registration
## Usage: script.sh <NP>
##   NP      participant id (e.g. sub-P001)
## Software: ANTs, FSL, AFNI
#################################################################################################################

set -uo pipefail

NP="$1"
tp=(TP1 TP2)
distmap_file="distancemap_core_eigenv_dwi.nii.gz"

sd_pre=/data/project/${tp[0]};
sd_post=/data/project/${tp[1]};
pref=${sd_pre}/derivatives/segm/${NP};
postf=${sd_post}/derivatives/segm/${NP};
outdir_reg=${sd_post}/derivatives/prepost_reg/${NP};
mkdir -p ${outdir_reg};

# inverse mask
for i in 0 1; do
fslmaths /data/project/${tp[i]}/derivatives/segm/${NP}/segm/segmentation.nii.gz -thr 1 -uthr 1 -bin mask1.nii.gz
fslmaths /data/project/${tp[i]}/derivatives/segm/${NP}/segm/segmentation.nii.gz -thr 3 -uthr 3 -bin mask3.nii.gz
fslmaths mask1.nii.gz -add mask3.nii.gz -bin /data/project/${tp[i]}/derivatives/segm/${NP}/segm/segmentation_core_bin.nii.gz
fslmaths /data/project/${tp[i]}/derivatives/segm/${NP}/segm/segmentation_core_bin.nii.gz -mul -1 -add 1 -mul /data/project/${tp[i]}/derivatives/segm/${NP}/preproc/t1c_bet_mask.nii.gz \
/data/project/${tp[i]}/derivatives/segm/${NP}/segm/mask_inv.nii.gz
rm mask1.nii.gz mask3.nii.gz
done

# register
t1_pre=${pref}/preproc/t1c_normalized_bet.nii.gz
mask_pre=${pref}/segm/mask_inv.nii.gz
t1_post=${postf}/preproc/t1c_normalized_bet.nii.gz
mask_post=${postf}/segm/mask_inv.nii.gz
antsRegistration \
--dimensionality 3 \
--float 0 \
--output [${outdir_reg}/post2pre_,${outdir_reg}/post2pre_Warped.nii.gz] \
--interpolation Linear \
--winsorize-image-intensities [0.005,0.995] \
--use-histogram-matching 1 \
--initial-moving-transform [${t1_pre},${t1_post},1] \
--transform Rigid[0.1] \
--metric MI[${t1_pre},${t1_post},1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] \
--shrink-factors 8x4x2x1 \
--smoothing-sigmas 4x2x1x0vox \
--transform Affine[0.1] \
--metric MI[${t1_pre},${t1_post},1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] \
--shrink-factors 8x4x2x1 \
--smoothing-sigmas 4x2x1x0vox \
--transform SyN[0.15,3,0] \
--metric CC[${t1_pre},${t1_post},1,4] \
--convergence [100x70x50x20,1e-6,10] \
--shrink-factors 8x4x2x1 \
--smoothing-sigmas 3x2x1x0vox \
-x [${mask_pre},${mask_post}]

# resample distance map
distmap=${sd_pre}/derivatives/DWI/${NP}/${distmap_file}
targ=${sd_post}/derivatives/GRE/${NP}/chisep/R2s.nii
mattarg=${sd_post}/derivatives/GRE/${NP}/reg/r2s2t1_0GenericAffine.mat
matt1cpost=${postf}/preproc/temp/atlas-space/M_atlas__t1c.mat
matpostpre=${sd_post}/derivatives/prepost_reg/${NP}/post2pre_0GenericAffine.mat
matt1c=${pref}/preproc/temp/atlas-space/M_atlas__t1c.mat
matdwi=${sd_pre}/derivatives/DWI/${NP}/reg/${NP}_dwi2t1_ants.mat

antsApplyTransforms \
-d 3 \
-i ${distmap} \
-r ${targ} \
-t [${mattarg}, 1 ] \
-t [${matt1cpost}, 1 ] \
-t [${matpostpre}, 1 ] \
-t ${matt1c} \
-t ${matdwi} \
-o ${outdir_reg}/${distmap_file/dwi.nii.gz/post_r2s.nii.gz} \
-n Linear

3dcalc -a ${outdir_reg}/${distmap_file/dwi.nii.gz/post_r2s.nii.gz} -expr 'ispositive(a)*1' -prefix binary_mask.nii.gz
3dmask_tool -input binary_mask.nii.gz -dilate_input -2 -prefix eroded_mask.nii.gz
3dcalc -a ${outdir_reg}/${distmap_file/dwi.nii.gz/post_r2s.nii.gz} -b eroded_mask.nii.gz -expr a*b -prefix ${outdir_reg}/${distmap_file/dwi.nii.gz/post_r2s_ero.nii.gz}
rm binary_mask.nii.gz eroded_mask.nii.gz

# resample tumor mask
antsApplyTransforms \
-d 3 \
-i ${mask_post} \
-r ${sd_post}/derivatives/GRE/${NP}/chisep/R2s.nii \
-t [${sd_post}/derivatives/GRE/${NP}/reg/r2s2t1_0GenericAffine.mat, 1 ] \
-t [${matt1cpost}, 1 ] \
-o ${postf}/segm/mask_inv_r2s.nii.gz \
-n NearestNeighbor
