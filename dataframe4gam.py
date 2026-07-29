import numpy as np
import pandas as pd
import nibabel as nib
import os

# --- Parameters ---
tp = "TP1"
BASE_DIR = f"/data/project/{tp}"
SCALAR_LIST = ["R2s", "ChiDia", "ChiPara", "ChiTot", "QSM"]  # add all scalar measures you want
dist_fname = "distancemap_core_eigenv_r2s_ero.nii.gz"

def process_participant_continuous(subj, measure, distance_fname):
    """
    Voxel-wise pipeline with global distance clipping for HFM maps.
    Returns DataFrame ready for GAM.
    """

    # --- Paths ---
    if tp == "TP1":
        subj_dir_dist = os.path.join(BASE_DIR, f"derivatives/DWI/{subj}")
    else:  # TP2
        subj_dir_dist = os.path.join(BASE_DIR, f"derivatives/prepost_reg/{subj}")
    subj_dir_dwi = os.path.join(BASE_DIR, f"derivatives/DWI/{subj}")
    subj_dir_gre = os.path.join(BASE_DIR, f"derivatives/GRE/{subj}")
    subj_dir_dose = os.path.join(BASE_DIR, f"derivatives/dose/{subj}")
    subj_dir_segm = os.path.join(BASE_DIR, f"derivatives/segm/{subj}")

    scalar_file = os.path.join(subj_dir_gre, f"chisep/{measure}.nii")  
    distance_file = os.path.join(subj_dir_dist, distance_fname)
    #dose_file = os.path.join(subj_dir_dose, "ProSoma_DICOM_RT_r2s_ero_rescaled.nii")
    tum_file_inv = os.path.join(subj_dir_segm, "mask_inv_r2s.nii.gz")
    vessel_dia_file = os.path.join(subj_dir_gre, "chisep/vesselMask_dia.nii")
    vessel_para_file = os.path.join(subj_dir_gre, "chisep/vesselMask_para.nii")
    mask_img = os.path.join(subj_dir_dwi,'reg/t1_fast_pve_2_r2s_thr_+edema.nii.gz')
    mask_dgm_img = os.path.join(BASE_DIR,f'derivatives/GRE/{subj}/reg/first/T1seg_all_fast_firstseg_bin_r2s.nii.gz')

    # --- Load maps ---
    scalar_map = nib.load(scalar_file).get_fdata()
    distance_map = nib.load(distance_file).get_fdata()
    #dose_map = nib.load(dose_file).get_fdata()
    tum_mask_inv = nib.load(tum_file_inv).get_fdata().astype(bool)
    vessel_dia_map = nib.load(vessel_dia_file).get_fdata().astype(bool)
    vessel_para_map = nib.load(vessel_para_file).get_fdata().astype(bool)
    mask_map = nib.load(mask_img).get_fdata().astype(bool)
    mask_dgm_map = nib.load(mask_dgm_img).get_fdata().astype(bool)

    # --- Valid voxels ---
    valid = (
        np.isfinite(scalar_map)
        & np.isfinite(distance_map)
        & (distance_map > 0)
        #& np.isfinite(dose_map)
        #& (dose_map > 0)
        & (vessel_dia_map == 0)
        & (vessel_para_map == 0)
        & tum_mask_inv
        & (mask_map == 1)
        & (mask_dgm_map == 0)
    )
    s_lo, s_hi = np.nanpercentile(scalar_map[valid], [0, 99])
    d_lo, d_hi = np.nanpercentile(distance_map[valid], [1, 95])

    valid_mask = (
        valid
        & (scalar_map > s_lo)
        & (scalar_map < s_hi)
        & (distance_map > d_lo)
        & (distance_map < d_hi)
    )

    # --- Export valid mask
    valid_mask_img = nib.Nifti1Image(valid_mask.astype(np.uint8), nib.load(scalar_file).affine)
    out_mask_file = os.path.join(subj_dir_gre, f"{subj}_valid_mask_r2s.nii.gz")
    nib.save(valid_mask_img, out_mask_file)
    print(f"Saved valid mask: {out_mask_file}")

    # Distance normalization
    valid_dist = distance_map[valid_mask]
    d_min = np.nanmin(valid_dist)
    d_max = np.nanmax(valid_dist)
    distance_norm_values = (valid_dist - d_min) / (d_max - d_min)
    distance_norm_values = np.clip(distance_norm_values, 0, 1)
    distance_norm = np.full(distance_map.shape, np.nan)
    distance_norm[valid_mask] = distance_norm_values

    # Dose normalization (same voxels)
    #valid_dose = dose_map[valid_mask]
    #dose_min = np.nanmin(valid_dose)
    #dose_max = np.nanmax(valid_dose)
    #dose_norm_values = (valid_dose - dose_min) / (dose_max - dose_min)
    #dose_norm_values = 1 - dose_norm_values  # inverted
    #dose_norm = np.full(dose_map.shape, np.nan)
    #dose_norm[valid_mask] = dose_norm_values

    # --- Normalize scalar ---
    scalar_norm = scalar_map.copy()
    mean_val = np.nanmean(scalar_map[valid_mask])
    std_val = np.nanstd(scalar_map[valid_mask])
    if std_val > 0:
        scalar_norm[valid_mask] = (scalar_map[valid_mask] - mean_val) / std_val

    # --- Extract voxel-wise data ---
    voxel_indices = np.where(valid_mask)

    df_voxel = pd.DataFrame({
        "participant": subj,
        "distance_raw": distance_map[voxel_indices],
        "distance_norm": distance_norm[voxel_indices],
        #"dose_raw": dose_map[voxel_indices],
        #"dose_norm": dose_norm[voxel_indices],
        "scalar_raw": scalar_map[voxel_indices],
        "scalar_norm": scalar_norm[voxel_indices],
        "max_distance": d_max,
        "min_distance": d_min
    })

    return df_voxel


# --- Main pipeline ---
if __name__ == "__main__":
    part_file = "/data/project/part.csv"
    df_part = pd.read_csv(part_file, sep="\t")
    participants = df_part['NP'].tolist()

    for measure in SCALAR_LIST:
        
        print(f"\nProcessing measure: {measure}")
        all_voxels = []

        for subj in participants:
            print(f"  Participant: {subj}")
            df_voxel = process_participant_continuous(
                subj=subj,
                measure=measure,
                distance_fname=dist_fname
            )
            all_voxels.append(df_voxel)

        df_voxels_out = pd.concat(all_voxels, ignore_index=True)

        # Save CSV
        out_csv = os.path.join(
            BASE_DIR,
            f"res/{measure}_distancemap_core_eigenv_r2s_ero_file.csv"
        )
        df_voxels_out.to_csv(out_csv, index=False)

        print(f"Saved CSV: {out_csv}")
