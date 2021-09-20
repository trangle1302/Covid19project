import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import io
import os
import json


SAMPLES_DEST = "/data/trang/covid19_data_CZ8746_annotation/segmentation"


if __name__ == "__main__":
    df_cells = pd.read_csv(os.path.join(SAMPLES_DEST, "cells_combined.csv")) 
    print(df_cells[['nuclei-size', 'cytosol-size', 'cell-size']].head())

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
    plt.savefig(os.path.join(SAMPLES_DEST, 'showme.png'))
    #plt.show()

    df_cells = df_cells[df_cells['cell-size']<7.5e4]
    df_cells = df_cells[df_cells['cell-size']>3500]
    df_cells.to_csv(os.path.join(SAMPLES_DEST, "cells_combined_nobigcell.csv"), index=False)

    print(df_cells['virus-cytosol-50to75percentmean'])


    