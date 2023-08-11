import os
import imageio
import numpy as np
from glob import glob
import pandas as pd

ACQUIRED_DATA_PATH = "/data/trang/HPA_CZ8780"#HPA_DV9903_Prescreen"
PLATEID = "CZ8780_plate_II"#"DV9903_240323_preHPA_II__2023-03-24T12_10_32-Measurement_1b"
SAMPLES_DEST = f"{ACQUIRED_DATA_PATH}_max_projection_annotation/segmentation/{PLATEID}"

def concat_df(samples):
    new_list = []
    cell_dfs_list = []
    i = 0
    for _, sample in enumerate(samples):
        try:
            cell_df = pd.read_csv(os.path.join(sample, 'cell_quantify.csv'))
            cell_dfs_list.append(cell_df)
        except Exception as e:
            new_list.append(sample)
            print(sample)
            i += 1
        cell_dfs_list.append(cell_df)
    if i > 0:
        print(f"Missing samples {new_list}")
    celld_fs_combined = pd.concat(cell_dfs_list)
    return celld_fs_combined
    

if __name__ == "__main__":
    samples = []
    for root, dirs, files in os.walk(SAMPLES_DEST):
        if ('cell_quantify.csv' in files):
            samples.append(root)
    samples.sort()
    celld_fs_combined = concat_df(samples)
    print(celld_fs_combined.head())
    print(celld_fs_combined.shape)
    save_path = os.path.join(SAMPLES_DEST,'cells_combined.csv')
    print(f"Saving to {save_path}")
    celld_fs_combined.to_csv(save_path, index=False)
