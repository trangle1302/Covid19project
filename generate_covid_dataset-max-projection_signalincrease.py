import os
import glob
import imageio
import configs as cfg

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
    SAVE_FOLDER = f"{cfg.IMG_FOLDER}_up"
    if not os.path.exists(SAVE_FOLDER):
        os.makedirs(SAVE_FOLDER)

    channels_ratio = {
        '_red.tiff': 5,
        '_green.tiff': 1,
        '_blue.tiff': 10,
        '_yellow.tiff': 5
    }
    intensity_increase(cfg.IMG_FOLDER, cfg.SAVE_FOLDER, channels_ratio)
