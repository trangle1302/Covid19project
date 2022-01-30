import os
import os
import glob
import imageio


def intensity_increase(image_folder, save_dir, channels_ratio):
    covid = glob.glob(image_folder + '/' + '*_red.tiff')
    imglist = [f[:-9] for f in covid]
    print(f'Found {len(imglist)} FOVs')

    # Multiply each channel
    for img_path in imglist:
        FOVname = os.path.basename(img_path)
        for ch in channels_ratio.keys():
            img = imageio.imread(img_path + ch)
            print(img_path + ch, img.max(), img.min())
            img[img == img.min()] = 0
            img = (img*channels_ratio[ch]).astype('uint16')

            print(os.path.join(save_dir, FOVname + ch), img.max(), img.min())
            imageio.imwrite(os.path.join(save_dir, FOVname + ch),img)


if __name__ == "__main__":
    # generate csvs, including all-tasks.csv, all-well-tasks.csv
    # where to host the csvs
    #IMG_FOLDER = "/data/trang/covid19_data_CZ8746_max_projection/10"
    #SEGMENTATION_FOLDER = "/data/trang/covid19_data_CZ8746_annotation/segmentation2"
    #CHANNELS = ["_blue.tiff", "_red.tiff", "_green.tiff", "_yellow.tiff"]
    IMG_FOLDER = "/data/trang/211111_COVID19_repurposing_Marianna_max_projection/11"
    SAVE_FOLDER = "/data/trang/211111_COVID19_repurposing_Marianna_max_projection/11_up"
    if not os.path.exists(SAVE_FOLDER):
        os.makedirs(SAVE_FOLDER)

    channels_ratio = {
        '_red.tiff': 5,
        '_green.tiff': 1,
        '_blue.tiff': 6,
        '_yellow.tiff': 5
    }
    intensity_increase(CSV_FOLDER, SAVE_FOLDER, channels_ratio)
