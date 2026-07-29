import os, sys
import numpy as np
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
sys.path.append('/data/users/q128ag/hfm_1/HFM_Python_Notebooks-master')
import FileIO  # Assuming you have the FileIO utility from the HFM packageHFM_Python_Notebooks-master 

def synthetizeMetric3(alpha, beta, Vx, Vy, Vz, dual=False):
    # Compute the squared norm of the eigenvector
    VNorm2 = Vx**2 + Vy**2 + Vz**2  # This is Vx^2 + Vy^2 + Vz^2
    
    # Check if the norm squared is zero (i.e., if any voxel has all components zero)
    # Avoid division by zero by setting c to zero in this case
    c = np.zeros_like(VNorm2, dtype=float)  # Default to zero (or np.nan if you prefer)
    mask_non_zero_norm = VNorm2 > 0  # Mask where norm squared is greater than 0

    # Compute the scaling factor c only where VNorm2 is non-zero
    e = 2 if dual else -2  # Eigenvalue relationship based on `dual`
    c[mask_non_zero_norm] = (alpha[mask_non_zero_norm]**e - beta[mask_non_zero_norm]**e) / VNorm2[mask_non_zero_norm]
    #c[mask_non_zero_norm] = (alpha**e - beta**e) / VNorm2[mask_non_zero_norm] #if alpha and beta are fixed

    # Compute the 6 components of the symmetric metric tensor
    metric = np.empty((Vx.shape[0], Vx.shape[1], Vx.shape[2], 6), dtype=float)

    # Compute tensor components using the safe `c`
    metric[:,:,:,0] = c * Vx**2 + beta  # Component 1 (Vx^2)
    metric[:,:,:,1] = c * Vx * Vy      # Component 2 (Vx*Vy)
    metric[:,:,:,2] = c * Vy**2 + beta # Component 3 (Vy^2)
    metric[:,:,:,3] = c * Vx * Vz      # Component 4 (Vx*Vz)
    metric[:,:,:,4] = c * Vy * Vz      # Component 5 (Vy*Vz)
    metric[:,:,:,5] = c * Vz**2 + beta # Component 6 (Vz^2)
    
    # Check SPD condition
    eps=1e-6
    eig1 = beta
    eig2 = beta + c * VNorm2
    min_eig = np.minimum(eig1, eig2)
    num_bad = np.sum(min_eig <= eps*0.5)
    if num_bad > 0:
        # warning: you might want to log or inspect the voxels where this happens
        print(f"Warning: {num_bad} voxels have non-positive minimal eigenvalue (<= {eps*0.5}).")
        # Optionally clamp final metric to a tiny positive definite fallback
        # For those voxels we can set metric = (eps, 0, eps, 0, 0, eps)  (diagonal eps)
        bad_mask = (min_eig <= eps*0.5)
        metric[bad_mask, :] = np.array([eps, 0.0, eps, 0.0, 0.0, eps])

    return metric

FileHFM_binary_dir = '~/hfm_1/HamiltonFastMarching-master/bin'
dir = '/data/project/TP1/derivatives'
sub = "sub-P001"

print("Processing:", sub)

# ---- Output file path ----
out_file = os.path.join(
    dir, f'DWI/{sub}/distancemap_core_eigenv_dwi.nii.gz'
)

# Load tum data
tum_img = nib.load(os.path.join(dir,f'segm/{sub}/segm/segmentation_dwi_core_bin.nii.gz'))
tum_data = tum_img.get_fdata()
seed_positions = np.array(np.where(tum_data > 0)).T

# Load mask data
mask_img = nib.load(os.path.join(dir,f'DWI/{sub}/reg/t1_fast_pve_2_dwi_thr_+edema.nii.gz'))
mask_data = mask_img.get_fdata()
mask_dgm_img = nib.load(os.path.join(dir,f'GRE/{sub}/reg/first/T1seg_all_fast_firstseg_bin_dwi.nii.gz'))
mask_dgm_data = mask_dgm_img.get_fdata()
mask_tot = mask_data.astype(bool) & ~mask_dgm_data.astype(bool)
mask_data_inverted = ~(mask_tot.astype(bool))

# Load temp data
temp_img = nib.load(os.path.join(dir,f'DWI/{sub}/reg/{sub}_b0_masked.nii'))
temp_data = temp_img.get_fdata()
temp_affine = temp_img.affine

# Load FOD data
fod_img = nib.load(os.path.join(dir,f'DWI/{sub}/{sub}_wmfod_3tissues_norm_peaks.nii.gz'))
fod_data = fod_img.get_fdata()
fod_data = np.nan_to_num(fod_data, nan=0, posinf=0, neginf=0)

# Load 1val data
epsilon = 1e-6
val_min = 0
val_max = 0.004
val_img = nib.load(os.path.join(dir, f'DWI/{sub}/{sub}_1val.nii.gz'))
val_data = val_img.get_fdata()
mask_data_bool = mask_data.astype(bool)
val_data = np.nan_to_num(val_data, nan=0, posinf=0, neginf=0)
val_data_norm = np.full_like(val_data, epsilon)
valid_voxels = mask_data_bool & (val_data > 0) & (val_data >= val_min) & (val_data <= val_max)
val_selected = val_data[valid_voxels]
val_selected_norm = (val_selected - val_min) / (val_max - val_min)
val_data_norm[valid_voxels] = val_selected_norm + epsilon

# Extract FODs
fod_max = np.nanmax(fod_data)
fod_data_normalized = fod_data / (fod_max + epsilon)
first_peak = fod_data_normalized[:, :, :, :3]
first_peak_magnitude = np.linalg.norm(first_peak, axis=-1)
second_peak = fod_data_normalized[:, :, :, 3:6]
second_peak_magnitude = np.linalg.norm(second_peak, axis=-1)
ratio_second_first = np.divide(second_peak_magnitude, first_peak_magnitude, where=(first_peak_magnitude > epsilon))

# Create eigenvalues
beta = np.log1p(ratio_second_first) / np.log1p(1) * 0.5 # Eigenvalue for the ortho direction
beta = np.clip(beta, 1e-8, None)
alpha = 1 - beta   # Eigenvalue for the primary plane
#alpha = 0.5
#beta = 0.5
alpha = np.nan_to_num(alpha, nan=0, posinf=0, neginf=0)
beta = np.nan_to_num(beta, nan=0, posinf=0, neginf=0)
#beta *= val_data_norm ##for the weighted
#alpha *= val_data_norm ##for the weighted

# Extract the eigenvector components from the 4D image
Vx = first_peak[:,:,:,0]  # First component (Vx) from the last dimension
Vy = first_peak[:,:,:,1]  # Second component (Vy) from the last dimension
Vz = first_peak[:,:,:,2]  # Third component (Vz) from the last dimension
Vnorm = np.sqrt(Vx**2 + Vy**2 + Vz**2) + 1e-8
Vx /= Vnorm
Vy /= Vnorm
Vz /= Vnorm

# Prepare HFM input
hfm_input = {}
hfm_input['arrayOrdering'] = 'RowMajor'  # Adjust to match the order of your data
hfm_input['dims'] = np.array(temp_data.shape)  # Set the dimensions of the FOD image
hfm_input['origin'] = np.array([0, 0, 0])  # Origin for your FOD data (adjust accordingly)
hfm_input['gridScale'] = 1.0  # Adjust grid scale if needed

# Set up seed and tip points for the Fast Marching solver
hfm_input['seeds'] = seed_positions  # Seeds: starting points for front propagation
hfm_input['walls'] = mask_data_inverted
hfm_input['exportValues'] = 1  # Whether to export computed values (distances)
hfm_input['sndOrder'] = 1      # Second order accuracy for the front propagation

# Riemannian metric, select if using "FileHFM_Riemann3" 
#hfm_input['dualMetric'] = synthetizeMetric3(alpha, beta, Vx, Vy, Vz, dual=True)

# Set speed if using "FileHFM_Isotropic3" : 1 = isotropic, val_data_norm = diffusion weighted
hfm_input['speed'] = val_data_norm

# Run the Fast Marching solver
#hfm_output = FileIO.WriteCallRead(hfm_input, "FileHFM_Riemann3", binary_dir=FileHFM_binary_dir)
hfm_output = FileIO.WriteCallRead(hfm_input, "FileHFM_Isotropic3", binary_dir=FileHFM_binary_dir)
distance_map=hfm_output["values"]
distance_map_nifti = nib.Nifti1Image(distance_map, affine=temp_affine)
nib.save(distance_map_nifti, out_file)