#!/bin/bash

#################################################################################################################
## Single patient / single timepoint intermodality registration: dwi -> t1 -> gre
## Usage: script.sh <tp> <NP>
##   tp      timepoint folder name under /data/project (e.g. TP1, TP2)
##   NP      participant id (e.g. sub-P001)
## Software: HD-BET, ANTs, FSL, MRtrix, Convert3D
#################################################################################################################

set -uo pipefail

tp="$1"
NP="$2"
sd=/data/project/${tp}; ##input main directory

### gre-t1 registration ###
anat_dir=${sd}/sourcedata/raws/${NP}/anat;
r2s_reg=${sd}/derivatives/GRE/${NP}/reg
mkdir -p ${r2s_reg};
t1c=${NP}_ce-gad_T1w.nii.gz
mov=${NP}_echo-1_part-mag_MEGRE.nii.gz

hd-bet -i ${anat_dir}/${mov} -o ${r2s_reg}/${mov/.nii.gz/_masked.nii.gz} -device cpu --save_bet_mask --disable_tta || true;
hd-bet -i ${anat_dir}/${t1c} -o ${r2s_reg}/${t1c/.nii.gz/_masked.nii.gz} -device cpu --save_bet_mask --disable_tta || true;

antsRegistration \
  --dimensionality 3 \
  --float 0 \
  --output [${r2s_reg}/r2s2t1_, ${r2s_reg}/r2s2t1_Warped.nii.gz] \
  --interpolation Linear \
  --winsorize-image-intensities [0.005,0.995] \
  --use-histogram-matching 1 \
  --initial-moving-transform [${r2s_reg}/${t1c/.nii.gz/_masked.nii.gz}, ${r2s_reg}/${mov/.nii.gz/_masked.nii.gz}, 1] \
  --transform Rigid[0.1] \
  --metric MI[${r2s_reg}/${t1c/.nii.gz/_masked.nii.gz}, ${r2s_reg}/${mov/.nii.gz/_masked.nii.gz}, 1, 32, Regular, 0.25] \
  --convergence [1000x500x250x100, 1e-6, 10] \
  --shrink-factors 8x4x2x1 \
  --smoothing-sigmas 3x2x1x0vox

### dwi-t1 registration ###
dwi_dir=${sd}/derivatives/DWI/${NP}
segm_dir=${sd}/derivatives/segm/${NP};
dwi_reg=${dwi_dir}/reg;
mkdir -p ${dwi_reg}
wm=${dwi_dir}/reg/t1_fast_pve_2.nii.gz
adc=${dwi_dir}/${NP}_adc.nii.gz
dwi_prep=${dwi_dir}/${NP}_upsampled.mif;
dwi_mask=${dwi_dir}/${NP}_mask_upsampled.mif;

dwiextract ${dwi_prep} - -bzero | mrmath - mean - -axis 3 | mrcalc - ${dwi_mask} -mult ${dwi_reg}/${NP}_b0_masked.nii;
mrconvert ${dwi_mask} ${dwi_mask/.mif/.nii.gz}
fast -t 1 -n 3 -H 0.1 -I 4 -l 20.0 -o ${dwi_reg}/t1_fast ${r2s_reg}/${t1c/.nii.gz/_masked.nii.gz}
fslmaths ${dwi_reg}/t1_fast_seg.nii.gz -thr 3 -uthr 3 -bin ${dwi_reg}/wm.nii.gz
epi_reg --epi=${dwi_reg}/${NP}_b0_masked.nii --t1=${t1c} --t1brain=${r2s_reg}/${t1c/.nii.gz/_masked.nii.gz} --out=${dwi_reg}/${NP}_dwi2t1 --wmseg=${wm}

### mask with edema in dwi/gre space ###
c3d_affine_tool -ref ${t1c} -src ${dwi_reg}/${NP}_b0_masked.nii ${dwi_reg}/${NP}_dwi2t1.mat -fsl2ras -oitk ${dwi_reg}/${NP}_dwi2t1_ants.mat
antsApplyTransforms -d 3 -i ${wm} -r ${dwi_reg}/${NP}_b0_masked.nii \
-t [${dwi_reg}/${NP}_dwi2t1_ants.mat, 1] -o ${wm/.nii.gz/_dwi.nii.gz} -n Linear
mrcalc ${adc} 0 -gt ${adc} 0.005 -lt -mul tmp_valid_mask.mif
mrcalc ${adc} tmp_valid_mask.mif -mul adc_clipped.mif
threshold=$(echo "$(mrthreshold adc_clipped.mif -mask ${wm/.nii.gz/_dwi.nii.gz} -percentile 50) + \
1.5 * ($(mrthreshold adc_clipped.mif -mask ${wm/.nii.gz/_dwi.nii.gz} -percentile 75) - \
$(mrthreshold adc_clipped.mif -mask ${wm/.nii.gz/_dwi.nii.gz} -percentile 25))" | bc -l)
mrthreshold adc_clipped.mif -abs $threshold -comparison le tmp_mask.mif
mrcalc tmp_mask.mif ${wm/.nii.gz/_dwi.nii.gz} -mul ${wm/.nii.gz/_dwi_thr.nii.gz}
fslmaths ${wm/.nii.gz/_dwi_thr.nii.gz} -thr 0.95 -bin ${wm/.nii.gz/_dwi_thr_bin.nii.gz}
rm -rf tmp_mask.mif tmp_valid_mask.mif adc_clipped.mif

antsApplyTransforms -d 3 -i ${segm_dir}/segm/segmentation.nii.gz -r ${dwi_reg}/${NP}_b0_masked.nii -t [${dwi_reg}/${NP}_dwi2t1_ants.mat, 1] -t [${segm_dir}/preproc/temp/atlas-space/M_atlas__t1c.mat, 1 ] -o ${segm_dir}/segm/segmentation_dwi.nii.gz -n GenericLabel
fslmaths ${segm_dir}/segm/segmentation_dwi.nii.gz -thr 2 -uthr 2 -bin ${segm_dir}/segm/segmentation_dwi_edema_bin.nii.gz
fslmaths ${segm_dir}/segm/segmentation_dwi.nii.gz -thr 1 -uthr 1 -bin tmp_1.nii.gz
fslmaths ${segm_dir}/segm/segmentation_dwi.nii.gz -thr 3 -uthr 3 -bin tmp_3.nii.gz
fslmaths tmp_1.nii.gz -add tmp_3.nii.gz ${segm_dir}/segm/segmentation_dwi_core_bin.nii.gz
fslmaths ${segm_dir}/segm/segmentation_dwi.nii.gz -bin ${segm_dir}/segm/segmentation_dwi_bin.nii.gz
rm tmp_1.nii.gz tmp_3.nii.gz
fslmaths ${wm/.nii.gz/_dwi_thr_bin.nii.gz} -add ${segm_dir}/segm/segmentation_dwi_edema_bin.nii.gz -bin ${wm/.nii.gz/_dwi_thr_+edema.nii.gz}
antsApplyTransforms -d 3 \
  -i ${wm/.nii.gz/_dwi_thr_+edema.nii.gz} \
  -r ${r2s_reg}/${mov/.nii.gz/_masked.nii.gz} \
  -t [${r2s_reg}/r2s2t1_0GenericAffine.mat,1] \
  -t ${dwi_reg}/${NP}_dwi2t1_ants.mat \
  -o [${wm/.nii.gz/_r2s_thr_+edema.nii.gz}, 0] \
  -n NearestNeighbor
antsApplyTransforms -d 3 \
  -i ${segm_dir}/segm/segmentation_dwi_edema_bin.nii.gz \
  -r ${r2s_reg}/${mov/.nii.gz/_masked.nii.gz} \
  -t [${r2s_reg}/r2s2t1_0GenericAffine.mat,1] \
  -t ${dwi_reg}/${NP}_dwi2t1_ants.mat \
  -o [${segm_dir}/segm/segmentation_r2s_edema_bin.nii.gz, 0] \
  -n NearestNeighbor

### deep grey matter segmentation in dwi/gre space ###
mkdir -p ${r2s_reg}/first
run_first_all -b -i ${r2s_reg}/${t1c/.nii.gz/_masked.nii.gz} -o ${r2s_reg}/first/T1seg
fslmaths ${r2s_reg}/first/T1seg_all_fast_firstseg.nii.gz -bin ${r2s_reg}/first/T1seg_all_fast_firstseg_bin.nii.gz
antsApplyTransforms -d 3 -i ${r2s_reg}/first/T1seg_all_fast_firstseg_bin.nii.gz -r ${dwi_reg}/${NP}_b0_masked.nii -t [${dwi_reg}/${NP}_dwi2t1_ants.mat, 1] -o ${r2s_reg}/first/T1seg_all_fast_firstseg_bin_dwi.nii.gz -n GenericLabel
antsApplyTransforms -d 3 -i ${r2s_reg}/first/T1seg_all_fast_firstseg_bin.nii.gz -r ${r2s_reg}/${mov/.nii.gz/_masked.nii.gz} -t [${r2s_reg}/r2s2t1_0GenericAffine.mat, 1] -o ${r2s_reg}/first/T1seg_all_fast_firstseg_bin_r2s.nii.gz -n GenericLabel

### dwi-gre transformation ###
matrices=("${r2s_reg}/r2s2t1_0GenericAffine.mat" "${dwi_reg}/${NP}_dwi2t1_ants.mat" )
for mat in "${matrices[@]}"; do
ConvertTransformFile 3 ${mat} ${mat/.mat/.txt}
transformconvert ${mat/.mat/.txt} itk_import ${mat/.mat/_mrtrix.txt}
done
transformcalc ${r2s_reg}/r2s2t1_0GenericAffine_mrtrix.txt invert ${r2s_reg}/t12r2s_0GenericAffine_mrtrix.txt
transformcompose ${dwi_reg}/${NP}_dwi2t1_ants_mrtrix.txt ${r2s_reg}/t12r2s_0GenericAffine_mrtrix.txt ${dwi_reg}/${NP}_dwi2r2s_mrtrix.txt
mrtransform ${dwi_dir}/${NP}_wmfod_3tissues_norm.mif ${dwi_dir}/${NP}_wmfod_3tissues_norm_r2s.mif -linear ${dwi_reg}/${NP}_dwi2r2s_mrtrix.txt -reorient_fod yes -template ${r2s_der}/chisep/R2s.nii
sh2peaks -nthreads ${n_threads} -threshold 0 ${dwi_dir}/${NP}_wmfod_3tissues_norm_r2s.mif ${dwi_dir}/${NP}_wmfod_3tissues_norm_r2s_peaks.nii.gz

antsApplyTransforms -d 3 \
  -i ${dwi_dir}/distancemap_core_eigenv_dwi.nii.gz \
  -r ${r2s_reg}/${mov/.nii.gz/_masked.nii.gz} \
  -t [${r2s_reg}/r2s2t1_0GenericAffine.mat,1] \
  -t ${dwi_reg}/${NP}_dwi2t1_ants.mat \
  -o [${dwi_dir}/distancemap_core_eigenv_r2s.nii.gz, 0] \
  -n Linear
3dcalc -a ${dwi_dir}/distancemap_core_eigenv_r2s.nii.gz -expr 'ispositive(a)*1' -prefix binary_mask.nii.gz
3dmask_tool -input binary_mask.nii.gz -dilate_input -2 -prefix eroded_mask.nii.gz
3dcalc -a ${dwi_dir}/distancemap_core_eigenv_r2s.nii.gz -b eroded_mask.nii.gz -expr a*b -prefix ${dwi_dir}/distancemap_core_eigenv_r2s_ero.nii.gz
rm -rf binary_mask.nii.gz eroded_mask.nii.gz
