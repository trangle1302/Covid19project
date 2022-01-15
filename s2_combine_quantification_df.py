import os
import imageio
import numpy as np
from glob import glob
import pandas as pd

SAMPLES_DEST = "/data/trang/211111_COVID19_repurposing_Marianna_max_projection_annotation/segmentation/12"

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
    celld_fs_combined.to_csv(os.path.join(SAMPLES_DEST,'cells_combined.csv'), index=False)
