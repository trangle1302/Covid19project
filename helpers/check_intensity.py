import os
import glob
from imageio import imread
from numpy.core.fromnumeric import mean
from numpy.lib.function_base import _place_dispatcher
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    CHANNELS = ["_blue.tiff", "_red.tiff", "_green.tiff", "_yellow.tiff"]
    PLATE_FOLDERS = "/data/hao/covid19_data_max_projection_combined"
    SAVE_FOLDERS = "/data/trang/211111_COVID19_repurposing_Marianna_max_projection"
    #PLATE_FOLDERS = "/data/trang/211111_COVID19_repurposing_Marianna_max_projection"
    plates = [p for p in os.listdir(PLATE_FOLDERS) if os.path.isdir(f"{PLATE_FOLDERS}/{p}")]
    
    if not os.path.exists(f'{SAVE_FOLDERS}/channel_intensity_plate1-8.csv'):
        intensity_df = []
        for plate in plates:
            imlist = glob.glob(f'{PLATE_FOLDERS}/{plate}/*red.tiff')
            imlist = [f[:-9] for f in imlist]
            for f in imlist:
                line = [f, str(plate)]
                try:
                    for ch in CHANNELS:
                        im = imread(f + ch)
                        print(f, im.max(), im.mean())
                        line += [str(im.max())]
                        line += [str(int(im.mean()))]
                    intensity_df += [line]
                except:
                    print(f'{f} lack some channel')
        column_names = ['filename','plate','blue_max', 'blue_mean','red_max', 'red_mean','green_max', 'green_mean','yellow_max', 'yellow_mean']
        intensity_df = pd.DataFrame(intensity_df, columns = column_names)
        intensity_df.to_csv(f'{SAVE_FOLDERS}/channel_intensity_plate1-8.csv', index=False)
    else:
        intensity_df = pd.read_csv(f'{SAVE_FOLDERS}/channel_intensity_plate1-8.csv')
        print(intensity_df.shape)
        print(intensity_df.plate.value_counts())

    for plate in plates:
        tmp = intensity_df[intensity_df.plate == str(plate)]
        print(f"Plate {plate}: {tmp.shape}")
        fig, ax = plt.subplots(2,4, figsize=(25,15))
        ax[0,0].hist(tmp.blue_max)
        ax[0,0].set_title('Max blue')
        ax[0,1].hist(tmp.red_max)
        ax[0,1].set_title('Max red')
        ax[0,2].hist(tmp.green_max)
        ax[0,2].set_title('Max green')
        ax[0,3].hist(tmp.yellow_max)
        ax[0,3].set_title('Max yellow')
        
        ax[1,0].hist(tmp.blue_max)
        ax[1,0].set_title('Mean blue')
        ax[1,1].hist(tmp.red_max)
        ax[1,1].set_title('Mean red')
        ax[1,2].hist(tmp.green_max)
        ax[1,2].set_title('Mean green')
        ax[1,3].hist(tmp.yellow_max)
        ax[1,3].set_title('Mean yellow')
        fig.savefig(f'{SAVE_FOLDERS}/Intensity_histogram_plate_{plate}.png')
