#!/bin/bash

#################################################################################################################
## Single patient / single timepoint DWI pipeline: preprocessing -> topup/eddy -> fods
## Usage: dwi_pipeline_single.sh <tp> <NP> <tum_hem>
##   tp      timepoint folder name under /data/project (e.g. TP1, TP2)
##   NP      participant id (e.g. sub-P001)
##   tum_hem tumour hemisphere, "lh" or "rh" (used for response-function masking in the fods step)
## Software: MRTrix, MRTrix3Tissue, SynB0, Freesurfer, FSL, AFNI
#################################################################################################################

set -uo pipefail

tp="$1"
NP="$2"
tum_hem="$3"

dir=/data/project/${tp}
path2freesurferlicense=~/FASTSURFER
path2synb0container=~/SYNB0
n_threads=16

healthy=$([[ "${tum_hem}" == "rh" ]] && echo "lh" || echo "rh")

dwi_raw=${dir}/sourcedata/raws/${NP}/dwi
dwi_der=${dir}/derivatives/DWI/${NP}
mkdir -p ${dwi_der}

data=${dwi_raw}/${NP}_dwi.nii.gz
bval=${data/nii.gz/bval}
bvec=${data/nii.gz/bvec}
json=${data/nii.gz/json}
t1=${dir}/sourcedata/raws/${NP}/anat/${NP}_T1w.nii.gz
pa_data=${dwi_raw}/${NP}_dwi_pa.nii.gz

cd "${dwi_der}" || exit 1

####################################################################################################################
############################################  PART 1 — PREPROCESSING  ##############################################
####################################################################################################################

echo -e "Participant ${NP} - preprocessing\n"

if [ ! -f "${NP}_dwi.mif" ]; then
mrconvert ${data} ${NP}_dwi.mif -fslgrad ${bvec} ${bval} -nthreads ${n_threads}
fi

### denoise ###
echo "########################################################################################"
echo -e "- denoise"
dwidenoise ${NP}_dwi.mif ${NP}_den.mif -noise ${NP}_noise.mif -nthreads ${n_threads}
mrcalc ${NP}_dwi.mif ${NP}_den.mif -subtract ${NP}_residuals.mif
echo "denoising done"
echo -e " -------> inspect ${dwi_der}/${NP}_noise.mif and ${dwi_der}/${NP}_residuals.mif\n"

### degibbs ###
echo "########################################################################################"
echo -e "- degibbs"
mrdegibbs ${NP}_den.mif ${NP}_degibbs.mif -nthreads ${n_threads} -axes 0,1 #plane slice acquisition: axial
echo -e "correction for gibbs artifact done\n"

### topup (synb0) ###
mkdir -p INPUTS
mkdir -p OUTPUTS

cp ${t1} INPUTS/T1.nii.gz

dwiextract ${NP}_degibbs.mif - -bzero | mrmath - mean INPUTS/b0.nii.gz -axis 3

temp1=$(awk -F" " '/TotalReadoutTime/{print $2}' ${json})
readout_time_dwi=$(echo "${temp1%,*}")

if [ -f "${pa_data}" ]; then
echo "0 -1 0 ${readout_time_dwi}" >> INPUTS/acqparams.txt
echo -e "0 1 0 0" >> INPUTS/acqparams.txt
topup --imain=INPUTS/b0.nii.gz --acqp=INPUTS/acqparams.txt --out=OUTPUTS/topup --fout=OUTPUTS/field_map.nii.gz --config=${FSLDIR}/etc/flirtsch/b02b0.cnf

else
echo "0 -1 0 ${readout_time_dwi}" >> INPUTS/acqparams.txt
echo -e "0 -1 0 0" >> INPUTS/acqparams.txt


echo "########################################################################################"
echo "run container for part ${NP}"
singularity run -e \
-B INPUTS/:/INPUTS \
-B OUTPUTS/:/OUTPUTS \
-B ${path2freesurferlicense}/license.txt:/extra/freesurfer/license.txt \
${path2synb0container}/synb0-disco_v3.0.sif

echo "SynB0 finished successfully"
fi

####################################################################################################################
###############################################  PART 2 — EDDY  ####################################################
####################################################################################################################

echo -e "Participant ${NP} - eddy\n"

### get mask ###
mrconvert ${NP}_degibbs.mif ${NP}_topup_eddy.nii -strides -1,+2,+3,+4 -export_grad_fsl eddy_bvecs eddy_bvals
applytopup --imain=${NP}_topup_eddy.nii --inindex=1 --datain=INPUTS/acqparams.txt --topup=OUTPUTS/topup --method=jac --out=${NP}_topupapplied
fslroi ${NP}_topupapplied.nii.gz b04mask.nii.gz 0 1
3dAutomask -prefix ${NP}_mask_posttopup.nii b04mask.nii.gz
rm -rf b04mask.nii.gz

### eddy ###
size=$(mrinfo -size ${NP}_degibbs.mif | awk '{print $4}')
for (( z=0; z<${size}; z++ )); do echo "1" >> ${NP}_index.txt; done

echo "########################################################################################"
echo "- eddy"
export OMP_NUM_THREADS=${n_threads}
eddy --imain=${NP}_topup_eddy.nii --mask=${NP}_mask_posttopup.nii --acqp=INPUTS/acqparams.txt --index=${NP}_index.txt \
--bvecs=eddy_bvecs --bvals=eddy_bvals --topup=OUTPUTS/topup --out=${NP}_post_eddy --repol

### post eddy ###
mrconvert ${NP}_post_eddy.nii.gz ${NP}_preproc.mif -strides -1,2,3,4 -fslgrad ${NP}_post_eddy.eddy_rotated_bvecs eddy_bvals
totalSlices=$(mrinfo ${NP}_preproc.mif | grep Dimensions | awk '{print $6 * $8}')
totalOutliers=$(awk '{ for(i=1;i<=NF;i++)sum+=$i } END { print sum }' ${NP}_post_eddy.eddy_outlier_map)
echo "If the following number is greater than 10, you may have to discard this subject because of too much motion or corrupted slices"
echo "scale=5; (${totalOutliers} / ${totalSlices} * 100)/1" | bc | tee ${NP}_percentageOutliers.txt
dwiextract ${NP}_preproc.mif - -bzero | mrmath - mean b04mask.nii -axis 3
3dAutomask -prefix ${NP}_posteddy_mask.nii b04mask.nii
rm -rf b04mask.nii ${NP}_mask_posttopup.nii ${NP}_topup_eddy.nii ${NP}_topupapplied.nii.gz ${NP}_topupapplied.mif eddy_bvals* eddy_bvecs* ${NP}_post_eddy.nii.gz ${NP}_index.txt

####################################################################################################################
###############################################  PART 3 — FODS  ####################################################
####################################################################################################################

### response function estimation (tumour: mask restricted to healthy hemisphere) ###
shells=$(tr ' ' '\n' < "${bval}" | sort -nu | wc -l)
if [ "${healthy}" == "rh" ]; then
3dcalc -a ${NP}_posteddy_mask.nii -dicom -expr 'isnegative(x-0)*a' -prefix ${NP}_mask_rh.nii.gz
elif [ "${healthy}" == "lh" ]; then
3dcalc -a ${NP}_posteddy_mask.nii -dicom -expr 'ispositive(x-0)*a' -prefix ${NP}_mask_lh.nii.gz
fi

echo "########################################################################################"
echo -e "- response function estimation"
dwi2response dhollander ${NP}_preproc.mif ${NP}_wmr.txt ${NP}_gmr.txt ${NP}_csfr.txt -voxels ${NP}_voxels_rf.mif -mask ${NP}_mask_${healthy}.nii.gz
echo -e "response function estimation done"
echo -e "check which voxels were used for estimation: cd ${dwi_der} && mrview ${NP}_preproc.mif -overlay.load ${NP}_voxels_rf.mif\n"

### upsampling ###
vox_size="1"

echo "########################################################################################"
echo -e "- upsampling"
mrgrid ${NP}_preproc.mif regrid ${NP}_upsampled.mif -voxel ${vox_size}
mrgrid ${NP}_posteddy_mask.nii regrid - -template ${NP}_upsampled.mif -interp linear -datatype bit | maskfilter - median ${NP}_mask_upsampled.mif
echo -e "upsampling done\n"

### fod calculation ###
outname="3tissues"
mrtrix3tissues="y"

echo "########################################################################################"
echo -e "- calculate fods"
if [ "${mrtrix3tissues}" == "y" ]; then
ss3t_csd_beta1 ${NP}_upsampled.mif ${NP}_wmr.txt ${NP}_wmfod_${outname}.mif ${NP}_gmr.txt ${NP}_gmfod_${outname}.mif \
${NP}_csfr.txt ${NP}_csffod_${outname}.mif -mask ${NP}_mask_upsampled.mif -nthreads ${n_threads}
echo "check fods: cd ${dwi_der} && fod2dec ${NP}_wmfod_${outname}.mif ${NP}_decfod.mif -mask ${NP}_mask_upsampled.mif && mrview ${NP}_decfod.mif -odf.load_sh ${NP}_wmfod_${outname}.mif"
else
if (( shells > 2 )); then
dwi2fod -mask ${NP}_mask_upsampled.mif msmt_csd ${NP}_upsampled.mif ${NP}_wmr.txt ${NP}_wmfod_${outname}.mif ${NP}_gmr.txt ${NP}_gmfod_${outname}.mif ${NP}_csfr.txt ${NP}_csffod_${outname}.mif
echo "check fods: cd ${dwi_der} && mrconvert -coord 3 0 ${NP}_wmfod_${outname}.mif - | mrcat ${NP}_csffod_${outname}.mif ${NP}_gmfod_${outname}.mif - ${NP}_vf.mif && mrview ${NP}_vf.mif -odf.load_sh ${NP}_wmfod_${outname}.mif"
else
dwi2fod -mask ${NP}_mask_upsampled.mif msmt_csd ${NP}_upsampled.mif ${NP}_wmr.txt ${NP}_wmfod_${outname}.mif ${NP}_csfr.txt ${NP}_csffod_${outname}.mif
echo "check fods: cd ${dwi_der} && mrconvert -coord 3 0 ${NP}_wmfod_${outname}.mif - | mrcat ${NP}_csffod_${outname}.mif - ${NP}_vf.mif && mrview ${NP}_vf.mif -odf.load_sh ${NP}_wmfod_${outname}.mif"
fi
fi
echo -e "fod calculation done\n"

### intensity normalization ###
echo "########################################################################################"
echo -e "- normalize fods"
mtnormalise ${NP}_wmfod_${outname}.mif ${NP}_wmfod_${outname}_norm.mif -mask ${NP}_mask_upsampled.mif
echo -e "global intensity normalization done\n"

### peaks extraction ###
echo "########################################################################################"
echo -e "- peaks extraction"
sh2peaks -nthreads ${n_threads} -threshold 0 ${NP}_wmfod_${outname}_norm.mif ${NP}_wmfod_${outname}_norm_peaks.nii.gz
echo -e "peaks extraction done\n"

### dti metrics ###
echo "########################################################################################"
echo -e "- calculate metrics"
dwi2tensor ${NP}_upsampled.mif -mask ${NP}_mask_upsampled.mif ${NP}_tensor.mif -nthreads ${n_threads}
tensor2metric ${NP}_tensor.mif -fa ${NP}_fa.nii.gz -adc ${NP}_adc.nii.gz -nthreads ${n_threads}
tensor2metric ${NP}_tensor.mif -value ${NP}_1val.nii.gz -nthreads ${n_threads}
echo -e "metrics calculated\n"


