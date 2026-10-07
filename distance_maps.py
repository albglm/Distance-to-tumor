"""
Geodesic distance-to-tumor maps with Hamiltonian Fast Marching (HFM)
=====================================================================

For every participant, computes the geodesic distance of each voxel of the
propagation domain from a seed region, by default the contrast-enhancing tumor
(CET), in diffusion space.

Distance maps (names = desc- label of the output file, used by
gam_distance_profiles.R):
    iso                 isotropic: uniform cost, distance in mm
    isoweighted         isotropic, cost scaled by axial diffusivity (primary map)
    aniso               anisotropic: cheaper along the main fiber direction
    anisoweighted       anisotropic and scaled by axial diffusivity
    contraisoweighted   isoweighted, seeded from the mirrored CET within the
                        contralateral hemisphere (contralateral control)
    isoweightedcetdil   isoweighted, seeded from the CET dilated by 2 voxels
                        (CET-boundary sensitivity analysis)
    isoweightedlesion   isoweighted, seeded from the whole lesion, i.e. distance
                        from the T2H-NAWM boundary (not used in the paper)
New maps are added as one line in MAPS (seeds, weighting, anisotropy). The
propagation model itself is in the block PROPAGATION MODEL (speed() and
anisotropy()), which can be edited to try other models.

Equations (see Methods):
    weighting   d = lambda1 / lambda_max + eps; voxels with lambda1 outside
                (0, lambda_max] get d = eps (effectively impassable)
    isotropic   |grad T| = 1 / F,  F = 1 (iso) or F = d (isoweighted)
    anisotropic grad T' D grad T = 1,  D = beta I + (alpha^2 - beta^2) u u'
                u = first FOD peak direction, r = second/first peak amplitude,
                beta = max(0.5 log2(1 + r), 1e-8), alpha = 1 - beta;
                alpha and beta multiplied by d in anisoweighted

HOW TO RUN
    1. Edit SETTINGS below (paths, maps to compute).
    2. python distance_maps.py                      all included participants
       python distance_maps.py sub-P001 sub-P002    only these participants
    Existing maps are skipped, so an interrupted run restarts where it stopped;
    set RECOMPUTE = True to overwrite them.

INPUT (BIDS-derivatives layout, diffusion space, 1 mm grid)
    <bids_root>/participants.tsv                       participant_id, include
    derivatives/masks/sub-<id>/dwi/
        sub-<id>_space-dwi_desc-tumor_dseg.nii.gz         tumor segmentation (TUMOR_LABELS);
                                                          seeds = CET + necrosis
        sub-<id>_space-dwi_label-NAWM_mask.nii.gz         thresholded FAST white matter
        sub-<id>_space-dwi_label-deepGM_mask.nii.gz       deep gray matter, excluded (optional)
        sub-<id>_space-dwi_label-CETmirrored_mask.nii.gz  contralateral control only
        sub-<id>_space-dwi_label-contrahemi_mask.nii.gz   contralateral control only
    derivatives/qmri/sub-<id>/dwi/
        sub-<id>_space-dwi_model-tensor_param-ad_dwimap.nii.gz   axial diffusivity (weighted maps only)
        sub-<id>_space-dwi_model-csd_param-peaks_dwimap.nii.gz   two largest FOD peaks, 6 volumes
                                                                 (anisotropic maps only)
    The output maps take the grid and orientation of the tumor segmentation.

OUTPUT
    derivatives/distancemaps/sub-<id>/dwi/sub-<id>_space-dwi_desc-<map>_distance.nii.gz

Requirements: Python 3 with numpy, nibabel, scipy; the HamiltonFastMarching
library (https://github.com/Mirebeau/HamiltonFastMarching), compiled locally.
Reference: Mirebeau J-M, Portegies J. Hamiltonian Fast Marching: a numerical
solver for anisotropic and non-holonomic eikonal PDEs. Image Processing On Line
2019;9:47-93. doi:10.5201/ipol.2019.227
"""

import csv
import os
import sys

import nibabel as nib
import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure

# ================================ SETTINGS ===================================

bids_root = "/path/to/bids"
hfm_python_dir = "/path/to/HamiltonFastMarching/Interfaces/PythonHFM/ExampleFiles/FileBased"
hfm_binary_dir = "/path/to/HamiltonFastMarching/bin"

MAPS_TO_RUN = ["iso", "isoweighted", "aniso", "anisoweighted"]
# contralateral control: ["contraisoweighted"]; CET sensitivity: ["isoweightedcetdil"]

# seeds: "cet" (CET + necrosis, default), "lesion" (CET + necrosis + T2H, i.e.
#        distance from the T2H-NAWM boundary) or "cetmirrored" (contralateral)
# seed_dilation: voxels (6-connected); contralateral: domain = contralateral hemisphere
MAPS = {
    "iso":               dict(anisotropic=False, weighted=False),
    "isoweighted":       dict(anisotropic=False, weighted=True),
    "aniso":             dict(anisotropic=True,  weighted=False),
    "anisoweighted":     dict(anisotropic=True,  weighted=True),
    "contraisoweighted": dict(anisotropic=False, weighted=True, seeds="cetmirrored", contralateral=True),
    "isoweightedcetdil": dict(anisotropic=False, weighted=True, seed_dilation=2),
    "isoweightedlesion": dict(anisotropic=False, weighted=True, seeds="lesion"),
}

# labels of the tumor segmentation (BraTS 2023 convention; BraTS 2021 uses cet=4)
TUMOR_LABELS = dict(necrosis=1, t2h=2, cet=3)

RECOMPUTE = False   # True = overwrite existing distance maps

LAMBDA_MAX = 4e-3   # mm^2/s, upper bound of the axial-diffusivity normalization
EPSILON = 1e-6

# ============================ PROPAGATION MODEL ==============================
# Edit these two functions to try another model; everything else stays the same.

def speed(ad, domain):
    """Weighting d = lambda1 / lambda_max + eps inside the domain; eps elsewhere
    and outside (0, lambda_max] (effectively impassable). Used by weighted maps."""
    d = np.full_like(ad, EPSILON)
    ok = domain & (ad > 0) & (ad <= LAMBDA_MAX)
    d[ok] = ad[ok] / LAMBDA_MAX + EPSILON
    return d


def anisotropy(r):
    """Speeds along (alpha) and across (beta) the first FOD peak, from the
    second-to-first peak amplitude ratio r (0 = single fiber, 1 = crossing).
    Used by anisotropic maps: D = beta I + (alpha^2 - beta^2) u u'."""
    beta = np.clip(0.5 * np.log2(1 + r), 1e-8, None)
    alpha = 1 - beta
    return alpha, beta

# =============================================================================

if not os.path.isdir(hfm_binary_dir):
    sys.exit(f"HFM binaries not found: {hfm_binary_dir} (set hfm_binary_dir in SETTINGS)")
sys.path.append(hfm_python_dir)
import FileIO  # noqa: E402  (HFM file interface)

# ================================ HELPERS ====================================

def deriv(pipeline, subject, name):
    return os.path.join(bids_root, "derivatives", pipeline, subject, "dwi",
                        f"{subject}_space-dwi_{name}.nii.gz")


def load(path):
    return np.nan_to_num(nib.load(path).get_fdata(), nan=0, posinf=0, neginf=0)


def fod_direction_and_ratio(path):
    """Unit vector of the first FOD peak and second-to-first peak amplitude ratio r."""
    peaks = load(path)
    first, second = peaks[..., 0:3], peaks[..., 3:6]
    first_amp = np.linalg.norm(first, axis=-1)
    second_amp = np.linalg.norm(second, axis=-1)
    u = first / (first_amp[..., None] + EPSILON)
    r = np.divide(second_amp, first_amp, out=np.zeros_like(first_amp), where=first_amp > EPSILON)
    return u, r


def dual_metric(alpha, beta, u):
    """D = beta I + (alpha^2 - beta^2) u u', stored as the 6 entries of a symmetric 3x3 tensor."""
    ux, uy, uz = u[..., 0], u[..., 1], u[..., 2]
    norm2 = ux**2 + uy**2 + uz**2
    c = np.zeros_like(norm2)
    ok = norm2 > 0
    c[ok] = (alpha[ok]**2 - beta[ok]**2) / norm2[ok]
    D = np.stack([c * ux**2 + beta, c * ux * uy, c * uy**2 + beta,
                  c * ux * uz, c * uy * uz, c * uz**2 + beta], axis=-1)
    # regularize voxels whose smallest eigenvalue is ~0
    bad = np.minimum(beta, beta + c * norm2) <= 0.5 * EPSILON
    D[bad] = [EPSILON, 0.0, EPSILON, 0.0, 0.0, EPSILON]
    return D

# ============================ ONE DISTANCE MAP ===============================

def output_file(subject, map_name):
    return os.path.join(bids_root, "derivatives", "distancemaps", subject, "dwi",
                        f"{subject}_space-dwi_desc-{map_name}_distance.nii.gz")


def distance_map(subject, map_name):
    spec = MAPS[map_name]
    mask = lambda label: load(deriv("masks", subject, f"label-{label}_mask")) > 0
    tumor_file = deriv("masks", subject, "desc-tumor_dseg")
    tumor = np.rint(load(tumor_file))
    cet = np.isin(tumor, [TUMOR_LABELS["necrosis"], TUMOR_LABELS["cet"]])   # CET + necrosis
    lesion = np.isin(tumor, list(TUMOR_LABELS.values()))                    # CET + necrosis + T2H

    # seeds
    seed_choice = spec.get("seeds", "cet")
    seeds = {"cet": lambda: cet, "lesion": lambda: lesion,
             "cetmirrored": lambda: mask("CETmirrored")}[seed_choice]()
    if spec.get("seed_dilation"):
        seeds = binary_dilation(seeds, generate_binary_structure(3, 1), iterations=spec["seed_dilation"])

    # propagation domain: NAWM + lesion, without deep gray matter
    domain = mask("NAWM") | lesion
    if os.path.exists(deriv("masks", subject, "label-deepGM_mask")):
        domain &= ~mask("deepGM")
    else:
        print(f"{subject}: no deep gray matter mask, nothing excluded")
    if spec.get("contralateral"):
        domain &= mask("contrahemi")

    # propagation cost
    d = speed(load(deriv("qmri", subject, "model-tensor_param-ad_dwimap")), domain) if spec["weighted"] else 1.0
    hfm_input = {
        "arrayOrdering": "RowMajor",
        "dims": np.array(domain.shape),
        "origin": np.array([0, 0, 0]),
        "gridScale": 1.0,
        "seeds": np.argwhere(seeds),
        "walls": ~domain,
        "exportValues": 1,
        "sndOrder": 1,
    }
    if spec["anisotropic"]:
        u, r = fod_direction_and_ratio(deriv("qmri", subject, "model-csd_param-peaks_dwimap"))
        alpha, beta = anisotropy(r)
        hfm_input["dualMetric"] = dual_metric(alpha * d, beta * d, u)
        solver = "FileHFM_Riemann3"
    else:
        hfm_input["speed"] = d
        solver = "FileHFM_Isotropic3"

    T = FileIO.WriteCallRead(hfm_input, solver, binary_dir=hfm_binary_dir)["values"]

    out = output_file(subject, map_name)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    nib.save(nib.Nifti1Image(T, nib.load(tumor_file).affine), out)
    print("saved:", out)

# ================================== MAIN =====================================

if __name__ == "__main__":
    subjects = sys.argv[1:]
    if not subjects:
        with open(os.path.join(bids_root, "participants.tsv")) as f:
            subjects = [row["participant_id"] for row in csv.DictReader(f, delimiter="\t")
                        if row["include"] == "1"]

    for subject in subjects:
        for map_name in MAPS_TO_RUN:
            if not RECOMPUTE and os.path.exists(output_file(subject, map_name)):
                print(f"{subject} {map_name}: exists, skipped")
                continue
            print(f"{subject} {map_name}: computing")
            try:
                distance_map(subject, map_name)
            except Exception as e:  # missing file, solver error: report and continue
                print(f"{subject} {map_name}: FAILED ({e})")
