"""
Synthetic example dataset
=========================

Writes a small BIDS-derivatives dataset of a few synthetic participants (a
spherical "brain" with a tumor in the left hemisphere) with every input file of
the pipeline, so the installation can be checked before using real data.

Each participant gets:
    - tumor segmentation (1 necrosis, 2 T2H, 3 CET), NAWM, hemisphere,
      deep gray matter and mirrored-CET masks
    - the six qMRI maps; their profiles change with distance from the CET and
      settle at the NAWM value after a small overshoot
    - axial diffusivity and FOD peaks (inputs of distance_maps.py)
    - Euclidean distance maps under all map names, so gam_distance_profiles.R
      can run without the HFM library

GRE and diffusion space are the same 1 mm grid here, so every file is written
in both spaces.

HOW TO RUN
    python make_example_data.py /path/to/example_bids
    then set bids_root to that folder in the other scripts. For a quick test
    set n_bootstrap <- 10 in gam_distance_profiles.R. To test distance_maps.py
    as well, set RECOMPUTE = True there (it overwrites the dwi distance maps).

Requirements: Python 3 with numpy, nibabel, scipy
"""

import os
import sys

import nibabel as nib
import numpy as np
from scipy.ndimage import distance_transform_edt

N_PARTICIPANTS = 3
SIZE = 48          # voxels per side, 1 mm
BRAIN_RADIUS = 21  # mm

# qMRI metrics: file name, NAWM value, profile amplitude, noise SD
METRICS = {
    "R2starmap":                     (20.0,  -8.0,  1.0),
    "Chimap":                        (-0.01, -0.01, 0.002),
    "desc-dia_Chimap":               (-0.02, -0.01, 0.002),
    "desc-para_Chimap":              (0.02,  -0.01, 0.002),
    "model-tensor_param-adc_dwimap": (8e-4,  1e-4,  2e-5),
    "model-tensor_param-fa_dwimap":  (0.45,  -0.25, 0.04),
}
MAP_NAMES = ["iso", "isoweighted", "aniso", "anisoweighted", "contraisoweighted"]


def save(root, pipeline, subject, space, name, data):
    folder = os.path.join(root, "derivatives", pipeline, subject, "dwi" if space == "dwi" else "anat")
    os.makedirs(folder, exist_ok=True)
    nib.save(nib.Nifti1Image(data.astype(np.float32), np.eye(4)),
             os.path.join(folder, f"{subject}_space-{space}_{name}.nii.gz"))


def participant(root, subject, rng):
    c = SIZE / 2
    x, y, z = np.meshgrid(*[np.arange(SIZE)] * 3, indexing="ij")
    brain = (x - c) ** 2 + (y - c) ** 2 + (z - c) ** 2 < BRAIN_RADIUS ** 2
    left = brain & (x < c)

    # tumor in the left hemisphere, irregular T2H
    t = np.array([c - 10, c, c]) + rng.integers(-2, 3, 3)
    r = np.sqrt((x - t[0]) ** 2 + (y - t[1]) ** 2 + (z - t[2]) ** 2)
    tumor = np.zeros(brain.shape)
    tumor[(r < 12 + 2 * np.sin(y / 3.0)) & brain] = 2
    tumor[r < 6] = 3
    tumor[r < 3] = 1
    cet = np.isin(tumor, [1, 3])
    mirrored = np.flip(cet, axis=0)

    deep_gm = (x - c - 5) ** 2 + (y - c) ** 2 + (z - c + 8) ** 2 < 16
    nawm = brain & (tumor == 0) & ~deep_gm & (rng.random(brain.shape) > 0.3)

    # profiles: deficit near the CET, small overshoot at about 5 mm, then NAWM value
    distance = distance_transform_edt(~cet)
    distance[~brain] = 0
    shape = np.exp(-distance / 1.5) - 0.15 * np.exp(-((distance - 5) / 2) ** 2)

    contra_distance = distance_transform_edt(~mirrored)
    contra_distance[~(brain & (x >= c))] = 0

    for space in ("gre", "dwi"):
        save(root, "masks", subject, space, "desc-tumor_dseg", tumor)
        save(root, "masks", subject, space, "label-NAWM_mask", nawm)
        save(root, "masks", subject, space, "label-tumorhemi_mask", left)
        save(root, "masks", subject, space, "label-deepGM_mask", deep_gm)
        save(root, "masks", subject, space, "label-CETmirrored_mask", mirrored)
        save(root, "masks", subject, space, "label-contrahemi_mask", brain & (x >= c))
        for name, (nawm_value, amplitude, noise) in METRICS.items():
            value = nawm_value + amplitude * shape + rng.normal(0, noise, brain.shape)
            save(root, "qmri", subject, space, name, np.where(brain, value, 0))
        for map_name in MAP_NAMES:
            d = contra_distance if map_name.startswith("contra") else distance
            save(root, "distancemaps", subject, space, f"desc-{map_name}_distance", d)

    peaks = np.zeros(brain.shape + (6,))
    peaks[..., 0], peaks[..., 4] = 1.0, 0.3            # first peak along x, second along y
    save(root, "qmri", subject, "dwi", "model-tensor_param-ad_dwimap", np.where(brain, 1.2e-3, 0))
    save(root, "qmri", subject, "dwi", "model-csd_param-peaks_dwimap", peaks)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python make_example_data.py /path/to/example_bids")
    root = sys.argv[1]
    os.makedirs(root, exist_ok=True)
    rng = np.random.default_rng(1)
    subjects = [f"sub-EX{i:02d}" for i in range(1, N_PARTICIPANTS + 1)]
    for s in subjects:
        participant(root, s, rng)
        print("written:", s)
    with open(os.path.join(root, "participants.tsv"), "w") as f:
        f.write("participant_id\tinclude\n" + "".join(f"{s}\t1\n" for s in subjects))
    print("example dataset in", root)
