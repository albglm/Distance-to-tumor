"""
GEODESIC DISTANCE MAP COMPUTATION, via the Hamiltonian Fast Marching method
==============================================================================
WHAT IT DOES
    For each patient, computes a "distance map": for every voxel (3D pixel)
    in a chosen brain region, how far it is from a chosen starting region. This
    isn't always straight-line distance -- how "far" a voxel counts as
    being can depend on the tissue it has to travel through (see map
    variants below). Used downstream by the GAM landmark analysis script.

MAP VARIANTS (set below)
    isotropic              - uniform cost everywhere -> true distance in mm
    isotropic, weighted     - cost scaled by axial diffusivity
    anisotropic             - cheaper to travel along the local fiber
                               direction than across it
    anisotropic, weighted   - anisotropic cost, also scaled by diffusivity

SEED AND REGION INPUTS ARE JUST MASKS
    Where the distance is measured from, and what tissue it's allowed to
    travel through, are both plain binary masks (see INPUT FILES below) --
    point them at whichever regions you want. E.g. `seed_mask_file` can be
    tumor core or tumor including edema; `allowed_region_mask_file` can
    be white matter near the tumor or a whole-brain mask; anything not in
    `allowed_region_mask_file` (or explicitly in `excluded_region_mask_file`,
    if used) is off-limits to the propagating front.

DEPENDENCIES
    Python: numpy, nibabel.
    External toolbox: HamiltonFastMarching, installed locally (https://github.com/Mirebeau/HamiltonFastMarching/tree/master).
 
REQUIRED INPUT FILES (per subject, NIfTI)
    A seed mask, an allowed-region mask, an optional excluded-region mask,
    a reference DWI (b=0) image, normalized FOD peaks, and a normalized
    axial-diffusivity ("1val") map, in the same space. See INPUT FILES below for details and
    which ones are only needed for certain MAP_TYPE/WEIGHTED settings.

OUTPUT
    <derivatives_dir>/DWI/<subject>/distancemap_core_<map_tag>_dwi.nii.gz
"""

import os
import sys
import numpy as np
import nibabel as nib  # reads/writes NIfTI (neuroimaging) files

sys.path.append('/toolboxes/HamiltonFastMarching/Interfaces/PythonHFM/ExampleFiles/FileBased')
import FileIO

# ============================== USER SETTINGS ================================

MAP_TYPE = "isotropic"   # "isotropic" or "anisotropic" -- see map variants above
WEIGHTED = False          # True: scale cost by axial diffusivity; False: unweighted

subject         = "sub-P001"
derivatives_dir  = "/data/project/TP1/derivatives"
hfm_binary_dir    = "/HamiltonFastMarching-master/bin"

VAL_MIN, VAL_MAX = 0, 0.004   # expected range of the axial-diffusivity map, for normalization
EPSILON = 1e-6                 # tiny constant to avoid divide-by-zero / exact-zero cost

# ================================ INPUT FILES =================================

seed_mask_file          = os.path.join(derivatives_dir, f"segm/{subject}/segm/tumor_core_dwi_bin.nii.gz")
    # Voxels > 0 here are the seeds the distance is measured from.
    # Point this at any binary mask instead -- e.g. tumor core or tumor including edema

allowed_region_mask_file = os.path.join(derivatives_dir, f"DWI/{subject}/reg/t1_fast_pve_2_dwi_thr_+edema.nii.gz")
    # Voxels > 0 here are the tissue the front is allowed to travel through.
    # Point this at any binary mask instead -- e.g. white matter mask or whole-brain mask

excluded_region_mask_file = os.path.join(derivatives_dir, f"GRE/{subject}/reg/t1_first_all_fast_firstseg_dwi_bin.nii.gz")
    # Optional: voxels > 0 here are removed from allowed_region_mask_file.
    # Set to None to skip this and use allowed_region_mask_file as-is.

reference_b0_file      = os.path.join(derivatives_dir, f"DWI/{subject}/{subject}_b0_masked.nii")  # for output grid size/orientation
fod_peaks_file         = os.path.join(derivatives_dir, f"DWI/{subject}/{subject}_wmfod_3tissues_norm_peaks.nii.gz")  # fiber orientation distribution (FOD) peaks; needed if MAP_TYPE = "anisotropic"
axial_diffusivity_file = os.path.join(derivatives_dir, f"DWI/{subject}/{subject}_1val.nii.gz")  # axial diffusivity (lambda1) map; needed if WEIGHTED = True

map_tag = MAP_TYPE + ("_weighted" if WEIGHTED else "")
out_file = os.path.join(derivatives_dir, f"DWI/{subject}/distancemap_core_{map_tag}_dwi.nii.gz")

def build_riemannian_metric(alpha, beta, Vx, Vy, Vz, dual=False):
    """
    Builds the per-voxel directional travel-cost tensor used for
    anisotropic propagation: an ellipsoid per voxel, cheap (`alpha`) along
    the local fiber direction (Vx,Vy,Vz), more costly (`beta`) across it.
    `dual=True` returns the convention expected by the solver's
    'dualMetric' input.
    """
    v_norm2 = Vx**2 + Vy**2 + Vz**2
    c = np.zeros_like(v_norm2, dtype=float)
    valid = v_norm2 > 0

    exponent = 2 if dual else -2
    c[valid] = (alpha[valid]**exponent - beta[valid]**exponent) / v_norm2[valid]

    metric = np.empty((*Vx.shape, 6), dtype=float)  # 6 = independent entries of a symmetric 3x3 tensor
    metric[..., 0] = c * Vx**2 + beta
    metric[..., 1] = c * Vx * Vy
    metric[..., 2] = c * Vy**2 + beta
    metric[..., 3] = c * Vx * Vz
    metric[..., 4] = c * Vy * Vz
    metric[..., 5] = c * Vz**2 + beta

    # A handful of voxels can end up with an invalid (non-positive) cost
    # ellipsoid due to numerical edge cases; replace those with a tiny,
    # harmless, direction-independent cost.
    eps = 1e-6
    min_eig = np.minimum(beta, beta + c * v_norm2)
    bad = min_eig <= eps * 0.5
    if bad.any():
        print(f"Warning: {bad.sum()} voxels have an invalid cost ellipsoid; using a fallback value.")
        metric[bad, :] = np.array([eps, 0.0, eps, 0.0, 0.0, eps])

    return metric


def load_normalized_diffusivity(path, mask, val_min=VAL_MIN, val_max=VAL_MAX, epsilon=EPSILON):
    """Loads the axial-diffusivity map and rescales in-mask, in-range voxels to (0, 1]."""
    val_data = np.nan_to_num(nib.load(path).get_fdata(), nan=0, posinf=0, neginf=0)
    val_norm = np.full_like(val_data, epsilon)
    valid = mask.astype(bool) & (val_data > 0) & (val_data >= val_min) & (val_data <= val_max)
    val_norm[valid] = (val_data[valid] - val_min) / (val_max - val_min) + epsilon
    return val_norm


def load_fod_direction(path, epsilon=EPSILON):
    """Loads the FOD peaks file and returns the dominant fiber direction per voxel, as a unit vector."""
    fod_data = np.nan_to_num(nib.load(path).get_fdata(), nan=0, posinf=0, neginf=0)
    Vx, Vy, Vz = fod_data[..., 0], fod_data[..., 1], fod_data[..., 2]
    v_norm = np.sqrt(Vx**2 + Vy**2 + Vz**2) + epsilon
    return Vx / v_norm, Vy / v_norm, Vz / v_norm


def get_fod_anisotropy_ratio(path, epsilon=EPSILON):
    """
    Loads the FOD peaks file and returns the second-peak / first-peak
    strength ratio per voxel: near 0 = one clear fiber direction, higher =
    crossing fibers / more isotropic. Sets how elongated each voxel's cost
    ellipsoid is (see build_riemannian_metric).
    """
    fod_data = np.nan_to_num(nib.load(path).get_fdata(), nan=0, posinf=0, neginf=0)
    fod_data = fod_data / (np.nanmax(fod_data) + epsilon)
    first_mag = np.linalg.norm(fod_data[:, :, :, :3], axis=-1)
    second_mag = np.linalg.norm(fod_data[:, :, :, 3:6], axis=-1)
    return np.divide(second_mag, first_mag, where=(first_mag > epsilon))


print("Processing:", subject, "| map type:", map_tag)

# --- seeds: where the distance is measured from -----------------------------
seed_data = nib.load(seed_mask_file).get_fdata()
seed_positions = np.array(np.where(seed_data > 0)).T

# --- walls: which tissue the front is allowed to travel through -------------
allowed_region = nib.load(allowed_region_mask_file).get_fdata().astype(bool)
if excluded_region_mask_file is not None:
    excluded_region = nib.load(excluded_region_mask_file).get_fdata().astype(bool)
    allowed_region = allowed_region & ~excluded_region
walls = ~allowed_region  # solver expects excluded voxels, not allowed ones

# --- reference geometry, so the output aligns with the subject's other scans ---
reference_img = nib.load(reference_b0_file)
grid_shape = reference_img.get_fdata().shape
reference_affine = reference_img.affine

# --- build the propagation cost: speed field (isotropic) or directional
#     cost tensor (anisotropic), scaled by diffusivity if WEIGHTED -------
if WEIGHTED:
    diffusivity_norm = load_normalized_diffusivity(axial_diffusivity_file, allowed_region)

if MAP_TYPE == "isotropic":
    speed = diffusivity_norm if WEIGHTED else 1  # 1 = uniform cost -> true geodesic distance in mm
    solver_name = "FileHFM_Isotropic3"

elif MAP_TYPE == "anisotropic":
    Vx, Vy, Vz = load_fod_direction(fod_peaks_file)
    anisotropy_ratio = get_fod_anisotropy_ratio(fod_peaks_file)

    # beta = cost across the fiber direction (grows with local isotropy); alpha = cost along it (kept low).
    beta = np.clip(np.log1p(anisotropy_ratio) / np.log1p(1) * 0.5, 1e-8, None)
    alpha = np.nan_to_num(1 - beta, nan=0, posinf=0, neginf=0)
    beta = np.nan_to_num(beta, nan=0, posinf=0, neginf=0)

    if WEIGHTED:
        alpha = alpha * diffusivity_norm
        beta = beta * diffusivity_norm

    dual_metric = build_riemannian_metric(alpha, beta, Vx, Vy, Vz, dual=True)
    solver_name = "FileHFM_Riemann3"

else:
    raise ValueError(f'MAP_TYPE must be "isotropic" or "anisotropic", got: {MAP_TYPE}')

hfm_input = {
    'arrayOrdering': 'RowMajor',
    'dims': np.array(grid_shape),
    'origin': np.array([0, 0, 0]),
    'gridScale': 1.0,
    'seeds': seed_positions,
    'walls': walls,
    'exportValues': 1,
    'sndOrder': 1,
}
if MAP_TYPE == "isotropic":
    hfm_input['speed'] = speed
else:
    hfm_input['dualMetric'] = dual_metric

hfm_output = FileIO.WriteCallRead(hfm_input, solver_name, binary_dir=hfm_binary_dir)

# =============================== SAVE OUTPUT ===================================

distance_map = hfm_output["values"]
nib.save(nib.Nifti1Image(distance_map, affine=reference_affine), out_file)
print("saved:", out_file)