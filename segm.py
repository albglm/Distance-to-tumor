from pathlib import Path
import shutil, os
from brainles_preprocessing.modality import Modality, CenterModality
from brainles_preprocessing.normalization.percentile_normalizer import (
    PercentileNormalizer,
)
from brainles_preprocessing.preprocessor import Preprocessor
from brainles_preprocessing.registration import ANTsRegistrator
from brainles_preprocessing.brain_extraction import HDBetExtractor
from brainles_preprocessing.defacing import QuickshearDefacer

import docker
client = docker.from_env()
from brats.core.segmentation_algorithms import AdultGliomaPreTreatmentSegmenter
from brats.core.segmentation_algorithms import AdultGliomaPreTreatmentAlgorithms
from brats import AdultGliomaPreAndPostTreatmentSegmenter
from brats.constants import AdultGliomaPreAndPostTreatmentAlgorithms

dir = Path("/data/project/TP1")
part_id = "sub-P001"
raw_data = dir / "sourcedata/raws" /part_id / "anat"
preproc_data = dir / "derivatives/segm" /part_id / "preproc"
preproc_data.mkdir(parents=True, exist_ok=True)
local_output = dir / "derivatives/segm" / part_id / "segm"
local_output.mkdir(parents=True, exist_ok=True)

t1c_file = raw_data / f"{part_id}_ce-gad_T1w.nii.gz"
t1_file = raw_data / f"{part_id}_T1w.nii.gz"
fla_file = raw_data / f"{part_id}_FLAIR.nii.gz"
t2_file = raw_data / f"{part_id}_T2w.nii.gz"


t1c_normalized_skull_output_path = preproc_data / "t1c_normalized_skull.nii.gz"
t1c_normalized_bet_output_path = preproc_data / "t1c_normalized_bet.nii.gz"
t1c_normalized_defaced_output_path = preproc_data / "t1c_normalized_defaced.nii.gz"
t1c_bet_mask = preproc_data / "t1c_bet_mask.nii.gz"
t1c_defacing_mask = preproc_data / "t1c_defacing_mask.nii.gz"

t1_normalized_bet_output_path = preproc_data / "t1_normalized_bet.nii.gz"
fla_normalized_bet_output_path = preproc_data / "fla_normalized_bet.nii.gz"
t2_normalized_bet_output_path = preproc_data / "t2_normalized_bet.nii.gz"

percentile_normalizer = PercentileNormalizer(
    lower_percentile=0.1,
    upper_percentile=99.9,
    lower_limit=0,
    upper_limit=1,
)

center = CenterModality(
    modality_name="t1c",
    input_path=t1c_file,
    normalizer=percentile_normalizer,
    # specify desired outputs, here we want to save the normalized skull, bet and defaced images
    normalized_skull_output_path=t1c_normalized_skull_output_path,
    normalized_bet_output_path=t1c_normalized_bet_output_path,
    normalized_defaced_output_path=t1c_normalized_defaced_output_path,
    # also save the bet and defacing mask
    bet_mask_output_path=t1c_bet_mask,
    defacing_mask_output_path=t1c_defacing_mask,
)

# Define the moving modalities, i.e. the modalities that are co-registered to the center modality.
moving_modalities = [
    Modality(
        modality_name="t1",
        input_path=t1_file,
        normalizer=percentile_normalizer,
        normalized_bet_output_path=t1_normalized_bet_output_path,
    ),
    Modality(
        modality_name="t2",
        input_path=t2_file,
        normalizer=percentile_normalizer,
        normalized_bet_output_path=t2_normalized_bet_output_path,
    ),
    Modality(
        modality_name="flair",
        input_path=fla_file,
        normalizer=percentile_normalizer,
        normalized_bet_output_path=fla_normalized_bet_output_path,
    ),
]

preprocessor = Preprocessor(
    center_modality=center,
    moving_modalities=moving_modalities,
    # Use ANTs for registration, other options are Niftyreg or eReg
    registrator=ANTsRegistrator(),
    # Use HDBet for brain extraction
    brain_extractor=HDBetExtractor(),
    # Use Quickshear for defacing,
    defacer=QuickshearDefacer(),
    # save temp
    temp_folder=preproc_data / "temp",
    use_gpu=True,
)

print(f"Preprocessing participant {part_id}...")
preprocessor.run()
print(f"Finished preprocessing {part_id}")

#segmentation
print(f"Segmenting participant {part_id}...")
segmenter = AdultGliomaPreAndPostTreatmentSegmenter(
    algorithm=AdultGliomaPreAndPostTreatmentAlgorithms.BraTS25_1,
    cuda_devices="0",
    force_cpu=False
)
    
segmenter.infer_single(
    t1c= preproc_data / "t1c_normalized_bet.nii.gz",
    t1n= preproc_data / "t1_normalized_bet.nii.gz",
    t2f= preproc_data / "fla_normalized_bet.nii.gz",
    t2w= preproc_data / "t2_normalized_bet.nii.gz",
    output_file= local_output / "segmentation.nii.gz"
)
print(f"Finished segmenting {part_id}")
