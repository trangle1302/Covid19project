import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import io
import os
import json
import configs as cfg

# Screening plate layout for Infected/Noninfected wells for different experiments
if cfg.EXPERIMENT == "CZ8751": # plate 11,12
    NONINFECTED = ["G11", "A12", "B12", "C12", "D12", "E12", "F12", "G12", "H12"]
    #INFECTED = 
elif cfg.EXPERIMENT == "CZ8746": # plate 10
    NONINFECTED = ["A","B","C","D"]
    INFECTED = ["E","F","G","H"]
elif cfg.EXPERIMENT == "DV9903":
    INFECTED = ["A","B","C","D","E","F","G"]
    NONINFECTED = ["H"]
elif cfg.EXPERIMENT == "CZ8780":
    NONINFECTED = ["B12", "C12", "D12", "E12", "H12"]


if __name__ == "__main__":
    df_cells = pd.read_csv(os.path.join(cfg.SAMPLES_DEST, "cells_combined.csv")) 
    print(df_cells[['nuclei-size', 'cytosol-size', 'cell-size']].head())
    print(f"Number of cells in total: {df_cells.shape[0]}")

    nuclei_sizes = df_cells['nuclei-size'].to_list()
    cytosol_sizes = df_cells['cytosol-size'].to_list()
    cell_sizes = df_cells['cell-size'].to_list()
    fig,axes = plt.subplots(nrows=1, ncols=3, figsize=(30,10))# , sharex=True, sharey=True
    ax = axes.ravel()
    ax[0].hist(nuclei_sizes, bins=500, range=(min(nuclei_sizes), 5E4))
    ax[0].set_title("nuclei_sizes", fontsize=40)
    ax[0].tick_params(axis='both', labelsize=20)
    ax[0].ticklabel_format(axis='both', style="sci", scilimits=(0,1))
    ax[1].hist(cytosol_sizes, bins=500, range=(min(cytosol_sizes), 2E5))
    ax[1].set_title("cytosol_sizes", fontsize=40)
    ax[1].tick_params(axis='both', labelsize=20)
    ax[1].ticklabel_format(axis='both', style="sci", scilimits=(0,1))
    ax[2].hist(cell_sizes, bins=500, range=(0, 1E5)) #min(cell_sizes)
    ax[2].set_title("cell_sizes", fontsize=40)
    ax[2].tick_params(axis='both', labelsize=20)
    ax[2].ticklabel_format(axis='both', style="sci", scilimits=(0,1))
    fig.tight_layout()
    plt.savefig(os.path.join(cfg.SAMPLES_DEST, 'showme.png'))
    #plt.show()
    print(df_cells["cell-size"].describe())
    df_cells = df_cells[df_cells['cell-size']<8e4]
    df_cells = df_cells[df_cells['cell-size']>3500]
    df_cells.to_csv(os.path.join(cfg.SAMPLES_DEST, "cells_combined_nobigcell.csv"), index=False)

    df_cells["paths"] = df_cells["image_id"]
    df_cells["image_id"] = [os.path.basename(f) for f in df_cells.image_id]
    df_cells["well_id"] = [os.path.basename(f).split("_")[-2] for f in df_cells.image_id]
    print(f"Well_id examples: {sorted(list(set(df_cells.well_id)))}")
    
    if cfg.EXPERIMENT == "CZ8751" or cfg.EXPERIMENT == "CZ8780": # plate 11,12
        df_cells["Infected"] = [0 if (f in NONINFECTED) else 1 for f in df_cells.well_id]
    else: # plate 10 and previous
        df_cells["Infected"] = [0 if (f[:1] in NONINFECTED) else 1 for f in df_cells.well_id]
    mean_noninfected = df_cells[df_cells["Infected"]==0]["virus-cytosol-50to75percentmean"].mean()
    df_cells["Infected"][(df_cells.Infected==1)&(df_cells["virus-cytosol-50to75percentmean"]<mean_noninfected)] = 0 
    print(f"NI/I threshold {mean_noninfected}, NI/I counts: {df_cells.Infected.value_counts()}")
    df_cells.to_csv(os.path.join(cfg.SAMPLES_DEST, "cells_combined_nobigcell.csv"), index=False)
    print(df_cells.groupby("well_id").agg({"Infected":"value_counts"}))
